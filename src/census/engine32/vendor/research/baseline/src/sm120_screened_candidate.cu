// Adapted from retained periodic candidate. Balanced gate added.
// CPU semantics tested; NO nvcc build or NVIDIA execution in this session.
#include "sm120_screened_api.cuh"
#include "two_cuts.hpp"
#include "../baseline/previous/src/warp_select.hpp"
#include "../baseline/previous/src/clipped_cover.hpp"
#include "../baseline/previous/baseline/vendor/src/wheel_filter.hpp"
namespace vdw_cut_gpu {
using namespace vdw_next;
constexpr unsigned FULL=0xffffffffu;
struct WarpPlan{vdw_periodic::DevicePlan p;u32 end;};
__device__ __forceinline__ WarpPlan build_warp_plan(vdw_periodic::Shape shape){
 u32 lane=threadIdx.x&31,P=shape.middle?vdw_periodic::period(shape.step):1;
 u32 sum=lane<P?(shape.middle?vdw_periodic::band_words(shape,lane):1):0;
 #pragma unroll
 for(u32 distance=1;distance<32;distance*=2){u32 before=__shfl_up_sync(FULL,sum,distance);if(lane>=distance)sum+=before;}
 u32 cycle=__shfl_sync(FULL,sum,P-1),tail=shape.middle%P;
 u32 tail_words=__shfl_sync(FULL,sum,tail?tail-1:0);if(!tail)tail_words=0;
 u32 mid=(shape.middle/P)*cycle+tail_words;
 return {{shape,P,cycle,mid,vdw_periodic::words(shape.begin)+mid+vdw_periodic::words(shape.end)},sum};
}
struct FetchPrefix{u32 value;__device__ __forceinline__ u32 operator()(u32 lane)const{return __shfl_sync(FULL,value,lane);}};
__device__ __forceinline__ vdw_periodic::Task select_pair(WarpPlan a,WarpPlan b,u32 rank){
 const unsigned mask=__ballot_sync(FULL,rank<a.p.total);
 // These branches are UNIFORM over the warp. Only a transition warp decodes twice.
 if(mask==FULL)return vdw_periodic::warp_select(a.p,rank,FetchPrefix{a.end});
 if(!mask)return vdw_periodic::warp_select(b.p,rank-a.p.total,FetchPrefix{b.end});
 auto x=vdw_periodic::warp_select(a.p,rank,FetchPrefix{a.end});
 auto y=vdw_periodic::warp_select(b.p,rank<a.p.total?b.p.total:rank-a.p.total,FetchPrefix{b.end});
 return rank<a.p.total?x:y;
}
template<int NP,int G>
__device__ __forceinline__ auto periodic_view(Input v,u32 row){
 if constexpr(G==1)return View{v.data,v.stride,v.B,NP,v.sign_xor};
 else{constexpr int NR=G==2?1:2;u32 c=row==0?1:G-1,mirror=NR-1-row;
  return ResidueView{v.data+row*NP*v.stride,v.data+mirror*NP*v.stride,v.stride,NP,v.sign_xor,(v.B-c)/G,(v.B-(G-c))/G};}
}
template<int NP,int G,bool SHARED>
__global__ void periodic_kernel(Input input,u32 t,u32 q,u32 D,const unsigned char*phase,
 SparseWord*out,unsigned long long capacity,unsigned long long*count,unsigned*overflow){
 constexpr int NR=G==1||G==2?1:2;constexpr u32 TAG=G==1?0:G==2?1:G==3?2:3;
 extern __shared__ u32 table[];Input local=input;local.np=NP;
 if constexpr(SHARED){for(u32 j=threadIdx.x;j<NR*NP*input.stride;j+=blockDim.x)table[j]=input.data[j];__syncthreads();local.data=table;}
 u32 lane=threadIdx.x&31,warp=(blockIdx.x*blockDim.x+threadIdx.x)/32,nwarps=gridDim.x*(blockDim.x/32);
 u32 slots=((D+5)/6)*(G==1||G==2?2:1)*NR,s=q-t+1;
 for(u32 job=warp;job<slots;job+=nwarps){u32 row=job%NR,n=job/NR,d;
  if constexpr(G==1)d=6*(n/2)+(n%2?5:1);if constexpr(G==2)d=6*(n/2)+(n%2?4:2);if constexpr(G==3)d=6*n+3;if constexpr(G==6)d=6*n+6;
  if(d>D)continue;u32 c=G==1?0:row==0?1:G-1,phi=phase?phase[d]:((t-1)/2)%s;
  if(phi>=s){if(!lane)atomicExch(overflow,1u);continue;}
  auto shapes=vdw_periodic::make_shapes(input.B,d,G,c,t,s,phi);
  auto a=build_warp_plan(shapes.a),b=build_warp_plan(shapes.b);u32 total=a.p.total+b.p.total;
  auto view=periodic_view<NP,G>(local,row);
  for(u32 base=0;base<total;base+=32){auto task=select_pair(a,b,base+lane);u32 bits=0;
   if(task.valid){
    bits=filter_range(view,d/G,t,task.word,task.first,task.last);
    // q is a fixed, witnessed FINE threshold for this epoch. Balanced two-cut
    // rejection keeps every anchor inside a coarse run of length >=q.
    bits=vdw_cut::screen(view,d/G,t,q,task.word,bits,1,nullptr);
   }
   unsigned mask=__ballot_sync(FULL,bits!=0);
   if(mask){unsigned long long start=0;if(!lane)start=atomicAdd(count,(unsigned long long)__popc(mask));start=__shfl_sync(FULL,start,0);
    if(bits){u32 rank=__popc(mask&((1u<<lane)-1u));auto i=start+rank;if(i<capacity)out[i]={d|(TAG<<30),G*32*task.word+c,bits,0};else atomicExch(overflow,1u);}
   }
  }
 }
}
template<int NP,int G,bool S>
cudaError_t do_launch(Input v,u32 t,u32 q,u32 D,const unsigned char*phase,SparseWord*out,unsigned long long capacity,unsigned long long*count,unsigned*overflow,int blocks,int threads,cudaStream_t stream){
 constexpr int NR=G==1||G==2?1:2;size_t bytes=S?(size_t)NR*NP*v.stride*4:0;
 if constexpr(S){int dev=0,maxbytes=0;auto e=cudaGetDevice(&dev);if(e!=cudaSuccess)return e;e=cudaDeviceGetAttribute(&maxbytes,cudaDevAttrMaxSharedMemoryPerBlockOptin,dev);if(e!=cudaSuccess)return e;if(bytes>(size_t)maxbytes)return cudaErrorInvalidValue;
  if(bytes>48*1024){e=cudaFuncSetAttribute(periodic_kernel<NP,G,S>,cudaFuncAttributeMaxDynamicSharedMemorySize,(int)bytes);if(e!=cudaSuccess)return e;}}
 periodic_kernel<NP,G,S><<<blocks,threads,bytes,stream>>>(v,t,q,D,phase,out,capacity,count,overflow);return cudaPeekAtLastError();
}
template<int NP,int G>
cudaError_t choice(Input v,u32 t,u32 q,u32 D,const unsigned char*phase,SparseWord*out,unsigned long long capacity,unsigned long long*count,unsigned*overflow,int blocks,int threads,bool shared,cudaStream_t stream){
 return shared?do_launch<NP,G,true>(v,t,q,D,phase,out,capacity,count,overflow,blocks,threads,stream):do_launch<NP,G,false>(v,t,q,D,phase,out,capacity,count,overflow,blocks,threads,stream);
}
cudaError_t launch_screened_group(Input v,u32 p,u32 g,u32 t,u32 q,u32 D,const unsigned char*phase,SparseWord*out,unsigned long long capacity,unsigned long long*count,unsigned*overflow,int blocks,int threads,bool shared,cudaStream_t stream){
 if(!v.data||!out||!count||!overflow||!capacity||blocks<1||blocks>65535||threads<32||threads>1024||threads%32||v.np<1||v.np>4||v.B<6||v.B>(1u<<20)||!(g==1||g==2||g==3||g==6)||t<3||q<t||q>256||!vdw_clipped::certified(p,{v.B,D,q})||v.stride<(v.B/g+32)/32+3||v.sign_xor>=(1u<<v.np))return cudaErrorInvalidValue;
 #define L(N,G) return choice<N,G>(v,t,q,D,phase,out,capacity,count,overflow,blocks,threads,shared,stream)
 #define GROUP(N) switch(g){case 1:L(N,1);case 2:L(N,2);case 3:L(N,3);case 6:L(N,6);}
 switch(v.np){case 1:GROUP(1);case 2:GROUP(2);case 3:GROUP(3);case 4:GROUP(4);}
 #undef GROUP
 #undef L
 return cudaErrorInvalidValue;
}
} // namespace
