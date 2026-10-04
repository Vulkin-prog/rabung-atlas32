#pragma once
#include <cuda_runtime.h>
#include "../vendor/src/sm120_candidate_api.cuh"
#include "anchor_device.hpp"
namespace vdw_next {
// Each g in {1,2,3,6} must be launched. No acceptance before all groups and checks.
// phase[d] is optional: null uses floor((t-1)/2) mod (q-t+1), also exact.
cudaError_t launch_anchor_group(Input input,u32 g,u32 t,u32 q,u32 D,
 const unsigned char* phase,SparseWord* output,unsigned long long capacity,
 unsigned long long* count,unsigned* overflow,int blocks,int threads,
 bool shared,cudaStream_t stream);
cudaError_t launch_anchor_completion(CompletionView v,const SparseWord* input,
 unsigned long long count,unsigned long long* best,unsigned* error,
 int blocks,int threads,cudaStream_t stream);
}
