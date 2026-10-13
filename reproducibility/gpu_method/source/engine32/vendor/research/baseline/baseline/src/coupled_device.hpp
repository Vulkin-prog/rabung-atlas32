#pragma once
// Exact local completion for a divisor profile and its refinement.
// CPU-compiled and tested. CUDA execution is NOT qualified by this header.
#include "../previous/baseline/vendor/src/completion.hpp"
#if defined(__CUDACC__)
#define VDW_PAIR_HD __host__ __device__ __forceinline__
#else
#define VDW_PAIR_HD inline
#endif
namespace vdw_pair {
using namespace vdw_next;
VDW_PAIR_HD int coarse_colour(int x,u32 coarse){
#if VDW_COLOUR_FAST
 // Current coupled profiles are 2->4 and 3->6; keep a general fallback.
 if(x<(int)(2*coarse))return x>=(int)coarse?x-(int)coarse:x;
#endif
 return x%(int)coarse;
}
using namespace vdw_next;
struct PairResult {u64 coarse_best,fine_best;u64 coarse_reads,fine_reads;bool ok;};
// Input colour entries are fine labels. Their reduction modulo coarse is exact.
// All marked anchors MUST be monochromatic for the coarse profile; this is the
// filter's contract. The precondition is tested separately against dense words.
VDW_PAIR_HD PairResult complete(CompletionView v,u32 coarse,SparseWord word,u64 cb,u64 fb){
 PairResult out{cb,fb,0,0,true};u32 d=word.tagged_d&0x3fffffffu,tag=word.tagged_d>>30;
 u32 g=tag==0?1:tag==1?2:tag==2?3:6;
 if(coarse<2||v.r>9||v.r<=coarse||v.r%coarse||v.t<3||v.t>=v.p||!v.colour||!v.inverse||word.right||!d||d>v.D||d>v.B||d%g){out.ok=false;return out;}
 u32 bits=word.left;
 while(bits){
#if defined(__CUDA_ARCH__)
  u32 k=__ffs(bits)-1;
#else
  u32 k=__builtin_ctz(bits);
#endif
  bits&=bits-1;u64 wide=(u64)word.first_b+(u64)g*k;
  if(!wide||wide>v.B){out.ok=false;return out;}u32 b=(u32)wide;
  if(!admissible(b,d,v.t))continue;
  u32 x=v.negative_quotient(b,d);
  if(!x||(u64)x+v.t>v.p){out.ok=false;return out;}
  i64 a=-(i64)b;int cf=v.quotient_at(a,x,d);if(cf<0){out.ok=false;return out;}int colour=coarse_colour(cf,coarse);
  u32 l=x,h=x+v.t-1;
  while(l>1){++out.coarse_reads;int c=v.quotient_at(a+((i64)l-1-x)*d,l-1,d);if(c<0){out.ok=false;return out;}if(coarse_colour(c,coarse)!=colour)break;--l;}
  while(h<v.p-1){++out.coarse_reads;int c=v.quotient_at(a+((i64)h+1-x)*d,h+1,d);if(c<0){out.ok=false;return out;}if(coarse_colour(c,coarse)!=colour)break;++h;}
  u32 reflected=v.p-h,start=l<reflected?l:reflected,len=h-l+1;
  u64 packed=pack_best(len,start);if(packed>out.coarse_best)out.coarse_best=packed;
  if(len<(out.fine_best>>32))continue;
  int previous=-1;u32 run=0;
  for(u32 pos=l;pos<=h;++pos){++out.fine_reads;int c=v.quotient_at(a+((i64)pos-x)*d,pos,d);if(c<0){out.ok=false;return out;}
   run=c==previous?run+1:1;previous=c;u32 st=pos-run+1,mirror=v.p-pos;if(mirror<st)st=mirror;
   packed=pack_best(run,st);if(packed>out.fine_best)out.fine_best=packed;
  }
 }
 return out;
}
} // namespace vdw_pair
#undef VDW_PAIR_HD
