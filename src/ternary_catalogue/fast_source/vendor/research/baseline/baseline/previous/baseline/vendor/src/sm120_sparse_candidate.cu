// Integration candidate. No GPU run or nvcc compilation is claimed in this delivery.
// Emits sparse WORD descriptors, not individual primitive fractions or final maxima.
#include <cuda_runtime.h>
#include "positive_filter.hpp"
#include "wheel_filter.hpp"
#include "sm120_candidate_api.cuh"
namespace vdw_next {
template<int NP,int G>
__device__ __forceinline__ auto make_view(Input v,u32 row){
 if constexpr(G==1)return View{v.data,v.stride,v.B,NP,v.sign_xor};
 else{constexpr int NR=G==2?1:2;u32 c=row==0?1:G-1,mirror=NR-1-row;
  return ResidueView{v.data+row*NP*v.stride,v.data+mirror*NP*v.stride,v.stride,NP,v.sign_xor,(v.B-c)/G,(v.B-(G-c))/G};}
}
template<int NP,int G,bool SHARED>
__global__ void emit_sparse(Input input,u32 t,u32 q,u32 D,SparseWord* output,
 unsigned long long capacity,unsigned long long* count,unsigned* overflow){
 constexpr int NR=G==1||G==2?1:2;
 extern __shared__ u32 scratch[];Input local=input;local.np=NP;
 if constexpr(SHARED){for(u32 i=threadIdx.x;i<NR*NP*input.stride;i+=blockDim.x)scratch[i]=input.data[i];__syncthreads();local.data=scratch;}
 const u32 lane=threadIdx.x&31;
 constexpr u32 MASK=0xffffffffu,GTAG=G==1?0:(G==2?1:(G==3?2:3));
 const u32 slots=((D+5)/6)*(G==1||G==2?2:1)*NR;
 for(u32 task=blockIdx.x;task<slots;task+=gridDim.x){
  u32 row=task%NR,n=task/NR,d=0,c=G==1?0:(row==0?1:G-1);
  if constexpr(G==1)d=6*(n/2)+(n%2?5:1);
  if constexpr(G==2)d=6*(n/2)+(n%2?4:2);
  if constexpr(G==3)d=6*n+3;
  if constexpr(G==6)d=6*n+6;
  if(d>D)continue;
  u32 bmin=((t-1)*d+1)/2,first=bmin>c?(bmin-c+G-1)/G:0,maximum=(input.B-c)/G;
  if(first>maximum)continue;
  auto view=make_view<NP,G>(local,row);u32 wfirst=first>>5,nw=(maximum>>5)-wfirst+1;
  for(u32 begin=0;begin<nw;begin+=blockDim.x){
   const u32 w=wfirst+begin+threadIdx.x;
   u32 bits=filter_range(view,d/G,t,w,first,maximum);
   Oriented out=prefilter_row(view,d/G,t,q,w,bits);
   const bool live=(out.left|out.right)!=0;
   // All 32 lanes participate. No gcd and no variable per-lane prefix in hot filter.
   const u32 mask=__ballot_sync(MASK,live);
   if(!mask)continue; // Warp-uniform empty case, no atomic or shuffle.
   unsigned long long base=0;
   if(lane==0)base=atomicAdd(count,(unsigned long long)__popc(mask));
   base=__shfl_sync(MASK,base,0);
   if(live){const u32 rank=__popc(mask&((1u<<lane)-1u));const unsigned long long index=base+rank;
    if(index<capacity)output[index]={d|(GTAG<<30),G*32*w+c,out.left,out.right};
    else atomicExch(overflow,1u);
   }
  }
 }
}
template<int NP,int G,bool SHARED>
cudaError_t launch_typed(Input v,u32 t,u32 q,u32 D,SparseWord*out,
 unsigned long long capacity,unsigned long long*count,unsigned*overflow,
 int blocks,int threads,cudaStream_t stream){
 constexpr int NR=G==1||G==2?1:2;
 size_t bytes=SHARED?(size_t)NR*NP*v.stride*sizeof(u32):0;
 if constexpr(SHARED){int device=0,maximum=0;auto e=cudaGetDevice(&device);if(e!=cudaSuccess)return e;
  e=cudaDeviceGetAttribute(&maximum,cudaDevAttrMaxSharedMemoryPerBlockOptin,device);if(e!=cudaSuccess)return e;
  if(bytes>(size_t)maximum)return cudaErrorInvalidValue;
  if(bytes>48*1024){e=cudaFuncSetAttribute(emit_sparse<NP,G,true>,cudaFuncAttributeMaxDynamicSharedMemorySize,(int)bytes);if(e!=cudaSuccess)return e;}
 }
 emit_sparse<NP,G,SHARED><<<blocks,threads,bytes,stream>>>(v,t,q,D,out,capacity,count,overflow);
 return cudaPeekAtLastError();
}
template<int NP,int G>
cudaError_t dispatch_shared(Input v,u32 t,u32 q,u32 D,SparseWord*out,
 unsigned long long capacity,unsigned long long*count,unsigned*overflow,
 int blocks,int threads,bool shared,cudaStream_t stream){
 return shared?launch_typed<NP,G,true>(v,t,q,D,out,capacity,count,overflow,blocks,threads,stream):launch_typed<NP,G,false>(v,t,q,D,out,capacity,count,overflow,blocks,threads,stream);
}
// Caller must reset count/overflow, launch EACH g=1,2,3,6 group, synchronize,
// reject overflow/error, and check all groups before any complete profile claim.
// Counts refer to SparseWord descriptors BEFORE primitive normalization.
cudaError_t launch_sparse_group(Input v,u32 g,u32 t,u32 q,u32 D,SparseWord*out,
 unsigned long long capacity,unsigned long long*count,unsigned*overflow,
 int blocks,int threads,bool shared,cudaStream_t stream){
 if(!v.data||!out||!count||!overflow||!capacity||blocks<1||threads<32||threads>1024||threads%32||v.np<1||v.np>4||v.B<6||v.B>(1u<<20)||!(g==1||g==2||g==3||g==6)||t<3||t>256||q<t||q>256||D==0||(unsigned long long)(t-1)*D>2ull*v.B||v.stride<(v.B/g+32)/32+3||v.sign_xor>=(1u<<v.np))return cudaErrorInvalidValue;
 #define GO(N,G) return dispatch_shared<N,G>(v,t,q,D,out,capacity,count,overflow,blocks,threads,shared,stream)
 #define GROUPS(N) switch(g){case 1:GO(N,1);case 2:GO(N,2);case 3:GO(N,3);case 6:GO(N,6);}
 switch(v.np){case 1:GROUPS(1);break;case 2:GROUPS(2);break;case 3:GROUPS(3);break;case 4:GROUPS(4);break;}
 #undef GROUPS
 #undef GO
 return cudaErrorInvalidValue;
}
} // namespace vdw_next
