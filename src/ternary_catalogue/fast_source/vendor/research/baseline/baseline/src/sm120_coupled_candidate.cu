// Candidate for the completion stage only. Not compiled or GPU-executed here.
#include <cuda_runtime.h>
#include "coupled_device.hpp"
namespace vdw_pair {
__global__ void complete_kernel(vdw_next::CompletionView v,unsigned coarse,
 const vdw_next::SparseWord* input,unsigned long long count,
 unsigned long long seed_coarse,unsigned long long seed_fine,
 unsigned long long* best_coarse,unsigned long long* best_fine,unsigned* error){
 unsigned long long i=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x;
 unsigned long long stride=(unsigned long long)gridDim.x*blockDim.x;
 for(;i<count;i+=stride){auto result=complete(v,coarse,input[i],seed_coarse,seed_fine);
  if(!result.ok){atomicExch(error,1u);continue;}
  atomicMax(best_coarse,result.coarse_best);atomicMax(best_fine,result.fine_best);
 }
}
// Preconditions not established by this wrapper: prime p, valid character table,
// inverses, every anchor monochromatic for coarse, and a complete generator epoch.
// Host orchestration must reject any upstream overflow, CUDA error, or error flag.
cudaError_t launch(vdw_next::CompletionView v,unsigned coarse,
 const vdw_next::SparseWord* input,unsigned long long count,
 unsigned long long seed_coarse,unsigned long long seed_fine,
 unsigned long long* best_coarse,unsigned long long* best_fine,unsigned* error,
 int blocks,int threads,cudaStream_t stream){
 if(!input||!best_coarse||!best_fine||best_coarse==best_fine||!error||!v.colour||!v.inverse||coarse<2||v.r<=coarse||v.r>9||v.r%coarse||blocks<1||threads<32||threads>1024||threads%32)
  return cudaErrorInvalidValue;
 if(!count)return cudaSuccess;
 complete_kernel<<<blocks,threads,0,stream>>>(v,coarse,input,count,seed_coarse,seed_fine,best_coarse,best_fine,error);
 return cudaPeekAtLastError();
}
}
