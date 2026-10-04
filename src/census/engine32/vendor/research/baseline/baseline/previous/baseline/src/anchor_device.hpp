#pragma once
#include "../vendor/src/completion.hpp"
#if defined(__CUDACC__)
#define VDW_AHD __host__ __device__ __forceinline__
#else
#define VDW_AHD inline
#endif
namespace vdw_next {
// Same 16-byte layout as SparseWord. left holds the anchor mask, right MUST be 0.
// This contract is NOT compatible with the old head-only completion.
VDW_AHD Completion complete_anchor_word(CompletionView v,SparseWord word,u64 seed){
 Completion out{seed,0,0,true};u32 d=word.tagged_d&0x3fffffffu,tag=word.tagged_d>>30;
 u32 g=tag==0?1:tag==1?2:tag==2?3:6;
 if(word.right||!d||d>v.D||d%g||d>v.B||!v.colour||!v.inverse){out.ok=false;return out;}
 u32 bits=word.left;
 while(bits){
#if defined(__CUDA_ARCH__)
  u32 k=__ffs(bits)-1;
#else
  u32 k=__builtin_ctz(bits);
#endif
  bits&=bits-1;u64 wide=(u64)word.first_b+(u64)g*k;
  if(!wide||wide>v.B){out.ok=false;return out;}
  u32 b=(u32)wide;if(!admissible(b,d,v.t))continue;++out.primitives;
  u32 x=(v.p-(u64)b*v.inverse[d]%v.p)%v.p;
  if(!x||(u64)x+v.t>v.p){out.ok=false;return out;}
  i64 a=-(i64)b;int col=v.small(a);if(col<0){out.ok=false;return out;}col=(col+v.r-v.colour[d])%v.r;
  u32 left=x,right=x+v.t-1;
  while(left>1){int c=v.quotient_at(a+((i64)left-1-x)*d,left-1,d);if(c<0){out.ok=false;return out;}if(c!=col)break;--left;}
  while(right<v.p-1){int c=v.quotient_at(a+((i64)right+1-x)*d,right+1,d);if(c<0){out.ok=false;return out;}if(c!=col)break;++right;}
  u32 reflected=v.p-right,start=left<reflected?left:reflected;
  u64 value=pack_best(right-left+1,start);if(value>out.best)out.best=value;
 }
 return out;
}
} // namespace vdw_next
#undef VDW_AHD
