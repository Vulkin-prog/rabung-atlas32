// Optional SECOND stage. Research candidate: neither nvcc-built nor GPU-tested
// in the recorded environment. Its input has already passed the coarse t-anchor.
#include <cuda_runtime.h>
#include "profile_masks.hpp"
namespace vdw_target_gpu {
using u32=unsigned int;
struct Row {
 vdw_next::View coarse_positive,fine_positive;
 vdw_next::ResidueView coarse_residue,fine_residue;
 u32 is_residue,selected_planes;
};
struct Input {u32 row,d,step,t,word,numerator_origin,coarse_anchor;};
struct Output {u32 d,numerator_origin,coarse_mask,fine_only_mask;};
// All rows and all input anchors in a launch belong to one threshold snapshot.
// Arrays of coarse/fine planes remain in global/L2 memory in this staged option;
// the original, expensive generator need not allocate additional shared memory.
// QC/QF are LOWER bounds backed by witnesses, never speculative estimates.
__global__ void gate(const Row*rows,u32 nr,const Input*input,u32 n,u32 qc,u32 qf,
                     Output*output,u32 capacity,u32*count,u32*failed){
 if((blockDim.x&31u)!=0){if(threadIdx.x==0)atomicExch(failed,1u);return;}
 const unsigned long long idx=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x;
 const u32 lane=threadIdx.x&31u;vdw_target::Masks result{0,0};Input task{};
 if(idx<n){task=input[idx];
  if(task.row>=nr||!task.step||task.t<3||qc<task.t||qf<task.t||qc<qf){atomicExch(failed,1u);}
  else{const Row row=rows[task.row];
   if(!row.selected_planes)atomicExch(failed,1u);
   else if(row.is_residue)result=vdw_target::split(row.coarse_residue,row.fine_residue,row.selected_planes,task.step,task.t,qc,qf,task.word,task.coarse_anchor);
   else result=vdw_target::split(row.coarse_positive,row.fine_positive,row.selected_planes,task.step,task.t,qc,qf,task.word,task.coarse_anchor);
  }
 }
 // No lane exits between entering this section and the matching collectives.
 const bool keep=(result.coarse|result.fine)!=0;
 const u32 votes=__ballot_sync(0xffffffffu,keep);
 if(votes){const u32 leader=__ffs(votes)-1u;u32 base=0;
  if(lane==leader)base=atomicAdd(count,(u32)__popc(votes));
  base=__shfl_sync(0xffffffffu,base,leader);
  if(keep){const u32 rank=__popc(votes&((1u<<lane)-1u));
   const unsigned long long dst=(unsigned long long)base+rank;
   if(dst>=capacity)atomicExch(failed,1u);
   else output[dst]={task.d,task.numerator_origin,result.coarse,result.fine};
  }
 }
}
} // namespace vdw_target_gpu
