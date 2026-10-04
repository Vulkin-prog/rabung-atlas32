// Device-compatible completion logic is CPU-tested; this CUDA wrapper is uncompiled.
#include <cuda_runtime.h>
#include "completion.hpp"
namespace vdw_next {
__global__ void finish_sparse(CompletionView v,const SparseWord* words,
 unsigned long long capacity,const unsigned long long* count,const unsigned* emission_overflow,
 unsigned long long seed,unsigned long long* best,unsigned* error){
 const unsigned long long i0=blockIdx.x*(unsigned long long)blockDim.x+threadIdx.x;
 const unsigned long long n=*count;
 if(*emission_overflow||n>capacity){if(i0==0)atomicExch(error,1u);return;}
 for(unsigned long long i=i0;i<n;i+=(unsigned long long)gridDim.x*blockDim.x){
  auto q=complete_word(v,words[i],seed);
  if(!q.ok){atomicExch(error,1u);continue;}
  if(q.best>seed)atomicMax(best,(unsigned long long)q.best);
 }
}
cudaError_t launch_sparse_completion(CompletionView v,const SparseWord*words,
 unsigned long long capacity,const unsigned long long*count,const unsigned*emission_overflow,
 unsigned long long seed,unsigned long long*best,unsigned*error,
 int blocks,int threads,cudaStream_t stream){
 u32 len=seed>>32,start=0xffffffffu-(u32)seed;
 if(!v.colour||!v.inverse||!words||!count||!emission_overflow||!best||!error||!capacity||blocks<1||threads<32||threads>1024||threads%32||v.p<7||!(v.p&1)||v.r<2||v.r>9||(v.p-1)%v.r||v.B<6||v.B>=v.p||v.B>(1u<<20)||v.t<3||v.t>256||v.t>=v.p||!v.D||(unsigned long long)(v.t-1)*v.D>2ull*v.B||!len||!start||(unsigned long long)start+len>v.p)return cudaErrorInvalidValue;
 // best MUST have been initialized to seed, and error to zero, before this launch.
 finish_sparse<<<blocks,threads,0,stream>>>(v,words,capacity,count,emission_overflow,seed,best,error);
 return cudaPeekAtLastError();
}
} // namespace vdw_next
