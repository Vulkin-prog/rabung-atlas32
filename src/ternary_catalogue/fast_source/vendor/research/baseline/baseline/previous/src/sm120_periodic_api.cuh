#pragma once
#include "../baseline/src/sm120_anchor_api.cuh"
namespace vdw_next {
// PRECONDITIONS: certified prime, certified Q-witness, correct input colour planes.
// All four groups must complete, no overflow/error, then bilateral completion.
// This wrapper checks the stronger clipped-volume sufficient condition.
cudaError_t launch_periodic_group(Input v,u32 p,u32 g,u32 t,u32 q,u32 D,
 const unsigned char*phase,SparseWord*out,unsigned long long capacity,
 unsigned long long*count,unsigned*overflow,int blocks,int threads,
 bool shared,cudaStream_t stream);
}
