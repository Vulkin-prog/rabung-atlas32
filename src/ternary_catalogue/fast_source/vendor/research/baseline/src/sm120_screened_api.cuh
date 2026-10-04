#pragma once
#include "../baseline/previous/src/sm120_periodic_api.cuh"
namespace vdw_cut_gpu {
// Required contracts: correct coarse planes, complete denominator epoch,
// prime p, a genuine fine q-witness, proper refinement, and clipped cover.
// The API cannot attest a witness supplied by a caller. Its host controller
// must verify that witness. All groups must complete, with zero error/overflow.
// Output is intentionally NOT a complete list of short anchors.
cudaError_t launch_screened_group(vdw_next::Input v,unsigned p,unsigned g,
 unsigned t,unsigned q,unsigned D,const unsigned char*phase,
 vdw_next::SparseWord*out,unsigned long long capacity,unsigned long long*count,
 unsigned*overflow,int blocks,int threads,bool shared,cudaStream_t stream);
}
