// Integration candidate, NOT compiled or executed on NVIDIA hardware in this session.
#include "sm120_anchor_api.cuh"
#include "../vendor/src/wheel_filter.hpp"
namespace vdw_next {
template<int NP,int G>
__device__ __forceinline__ auto anchor_view(Input v,u32 row){
 if constexpr(G==1)return View{v.data,v.stride,v.B,NP,v.sign_xor};
 else{constexpr int NR=G==2?1:2;u32 c=row==0?1:G-1,mirror=NR-1-row;
  return ResidueView{v.data+row*NP*v.stride,v.data+mirror*NP*v.stride,v.stride,NP,v.sign_xor,(v.B-c)/G,(v.B-(G-c))/G};}
}
template<int G,class V>
__device__ __forceinline__ void anchor_band(V view,u32 d,u32 t,u32 c,u32 blo,u32 bhi,
 SparseWord*output,unsigned long long capacity,unsigned long long*count,unsigned*overflow){
 constexpr unsigned FULL=0xffffffffu;constexpr u32 TAG=G==1?0:G==2?1:G==3?2:3;
 if(blo>bhi||bhi<c)return;u32 first=blo>c?(blo-c+G-1)/G:0,last=(bhi-c)/G;
 if(first>last)return;u32 lane=threadIdx.x&31,w0=first>>5,nw=(last>>5)-w0+1;
 // Every lane executes the same loop bounds and every warp collective.
 for(u32 chunk=0;chunk<nw;chunk+=32){u32 w=w0+chunk+lane;
  u32 bits=filter_range(view,d/G,t,w,first,last);unsigned mask=__ballot_sync(FULL,bits!=0);
  if(mask){unsigned long long base=0;if(!lane)base=atomicAdd(count,(unsigned long long)__popc(mask));base=__shfl_sync(FULL,base,0);
   if(bits){unsigned rank=__popc(mask&((1u<<lane)-1u));auto index=base+rank;
    if(index<capacity)output[index]={d|(TAG<<30),G*32*w+c,bits,0};else atomicExch(overflow,1u);}
  }
 }
}
template<int NP,int G,bool SHARED>
__global__ void emit_anchors(Input input,u32 t,u32 q,u32 D,const unsigned char* phases,
 SparseWord*output,unsigned long long capacity,unsigned long long*count,unsigned*overflow){
 constexpr int NR=G==1||G==2?1:2;extern __shared__ u32 scratch[];Input local=input;local.np=NP;
 if constexpr(SHARED){for(u32 j=threadIdx.x;j<NR*NP*input.stride;j+=blockDim.x)scratch[j]=input.data[j];__syncthreads();local.data=scratch;}
 u32 warp=(blockIdx.x*blockDim.x+threadIdx.x)/32,nwarps=gridDim.x*(blockDim.x/32);
 u32 slots=((D+5)/6)*(G==1||G==2?2:1)*NR,s=q-t+1;
 for(u32 task=warp;task<slots;task+=nwarps){u32 row=task%NR,n=task/NR,d;
  if constexpr(G==1)d=6*(n/2)+(n%2?5:1);if constexpr(G==2)d=6*(n/2)+(n%2?4:2);
  if constexpr(G==3)d=6*n+3;if constexpr(G==6)d=6*n+6;
  if(d>D)continue;u32 c=G==1?0:row==0?1:G-1,minimum=((t-1)*d+1)/2;
  auto view=anchor_view<NP,G>(local,row);
  if(d<=2){anchor_band<G>(view,d,t,c,minimum,input.B,output,capacity,count,overflow);continue;}
  u32 phase=phases?phases[d]:((t-1)/2)%s;
  if(phase>=s){if((threadIdx.x&31)==0)atomicExch(overflow,1u);continue;}
  u32 other=(t-2+s-phase)%s,first_n=minimum/d,max_n=input.B/d;
  // Jump DIRECTLY between selected bands; do not traverse omitted numerator words.
  for(u32 side=0;side<(phase==other?1u:2u);++side){u32 phi=side?other:phase;
   u32 n0=first_n+(phi+s-first_n%s)%s;
   for(u32 k=n0;k<=max_n;k+=s){u32 base=k*d;
    u32 lo=base+(side?d/2+1:1),hi=base+(phase==other?d-1:side?d-1:(d-1)/2);
    lo=lo<minimum?minimum:lo;hi=hi>input.B?input.B:hi;
    anchor_band<G>(view,d,t,c,lo,hi,output,capacity,count,overflow);
   }
  }
 }
}
template<int NP,int G,bool SHARED>
cudaError_t anchor_launch(Input input,u32 t,u32 q,u32 D,const unsigned char*phase,SparseWord*out,
 unsigned long long capacity,unsigned long long*count,unsigned*overflow,int blocks,int threads,cudaStream_t st){
 constexpr int NR=G==1||G==2?1:2;size_t bytes=SHARED?(size_t)NP*NR*input.stride*4:0;
 if constexpr(SHARED){int dev=0,maxbytes=0;auto e=cudaGetDevice(&dev);if(e!=cudaSuccess)return e;e=cudaDeviceGetAttribute(&maxbytes,cudaDevAttrMaxSharedMemoryPerBlockOptin,dev);if(e!=cudaSuccess)return e;if(bytes>(size_t)maxbytes)return cudaErrorInvalidValue;
 if(bytes>48*1024){e=cudaFuncSetAttribute(emit_anchors<NP,G,true>,cudaFuncAttributeMaxDynamicSharedMemorySize,(int)bytes);if(e!=cudaSuccess)return e;}}
 emit_anchors<NP,G,SHARED><<<blocks,threads,bytes,st>>>(input,t,q,D,phase,out,capacity,count,overflow);return cudaPeekAtLastError();
}
template<int NP,int G>
cudaError_t anchor_select(Input v,u32 t,u32 q,u32 D,const unsigned char*phase,SparseWord*out,unsigned long long capacity,unsigned long long*count,unsigned*overflow,int blocks,int threads,bool shared,cudaStream_t st){
 return shared?anchor_launch<NP,G,true>(v,t,q,D,phase,out,capacity,count,overflow,blocks,threads,st):anchor_launch<NP,G,false>(v,t,q,D,phase,out,capacity,count,overflow,blocks,threads,st);
}
cudaError_t launch_anchor_group(Input v,u32 g,u32 t,u32 q,u32 D,const unsigned char*phase,SparseWord*out,unsigned long long capacity,unsigned long long*count,unsigned*overflow,int blocks,int threads,bool shared,cudaStream_t st){
 if(!v.data||!out||!count||!overflow||!capacity||blocks<1||blocks>65535||threads<32||threads>1024||threads%32||v.np<1||v.np>4||v.B<6||v.B>(1u<<20)||!(g==1||g==2||g==3||g==6)||t<3||q<t||q>256||!D||(unsigned long long)(q-1)*D>2ull*v.B||v.stride<(v.B/g+32)/32+3||v.sign_xor>=(1u<<v.np))return cudaErrorInvalidValue;
 #define LAUNCH(N,G) return anchor_select<N,G>(v,t,q,D,phase,out,capacity,count,overflow,blocks,threads,shared,st)
 #define GROUP(N) switch(g){case 1:LAUNCH(N,1);case 2:LAUNCH(N,2);case 3:LAUNCH(N,3);case 6:LAUNCH(N,6);}
 switch(v.np){case 1:GROUP(1);case 2:GROUP(2);case 3:GROUP(3);case 4:GROUP(4);}
 #undef GROUP
 #undef LAUNCH
 return cudaErrorInvalidValue;
}
__global__ void complete_anchors(CompletionView view,const SparseWord*data,unsigned long long count,unsigned long long*best,unsigned*error){
 u64 local=0;bool ok=true;
 for(unsigned long long i=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x;i<count;i+=(unsigned long long)blockDim.x*gridDim.x){auto result=complete_anchor_word(view,data[i],local);if(!result.ok){ok=false;break;}local=result.best;}
 if(!ok)atomicExch(error,1u);if(local)atomicMax(best,(unsigned long long)local);
}
cudaError_t launch_anchor_completion(CompletionView v,const SparseWord*input,unsigned long long count,unsigned long long*best,unsigned*error,int blocks,int threads,cudaStream_t st){
 if(!input||!best||!error||!v.colour||!v.inverse||!v.D||v.B>=v.p||v.t<3||v.r<2||v.r>9||(v.p-1)%v.r||blocks<1||blocks>65535||threads<32||threads>1024||threads%32)return cudaErrorInvalidValue;
 complete_anchors<<<blocks,threads,0,st>>>(v,input,count,best,error);return cudaPeekAtLastError();
}
} // namespace vdw_next
