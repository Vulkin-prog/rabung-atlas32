#pragma once
#include <cuda_runtime.h>
#include "sparse_words.hpp"
namespace vdw_next {
struct Input {const std::uint32_t* data;std::uint32_t stride,B,np,sign_xor;};
cudaError_t launch_sparse_group(Input view,std::uint32_t g,std::uint32_t t,
 std::uint32_t q,std::uint32_t D,SparseWord* output,unsigned long long capacity,
 unsigned long long* count,unsigned* overflow,int blocks,int threads,
 bool shared,cudaStream_t stream);
}
