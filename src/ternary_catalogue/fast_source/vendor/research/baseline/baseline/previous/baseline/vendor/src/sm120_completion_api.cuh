#pragma once
#include <cuda_runtime.h>
#include "completion.hpp"
namespace vdw_next {
cudaError_t launch_sparse_completion(CompletionView view,const SparseWord* words,
 unsigned long long capacity,const unsigned long long* count,const unsigned* emission_overflow,
 unsigned long long seed,unsigned long long* best,unsigned* error,
 int blocks,int threads,cudaStream_t stream);
}
