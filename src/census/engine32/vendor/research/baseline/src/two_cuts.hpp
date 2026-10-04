#pragma once
// Certified rejection of a short COARSE run containing an interior anchor.
// The input mask has passed the t-symbol anchor filter. Unknown entries stay
// admissible. There is no head-only assumption. No extra shared-memory array.
#include "../baseline/previous/baseline/vendor/src/wheel_filter.hpp"
#if defined(__CUDACC__)
#define VDW_CUT_HD __host__ __device__ __forceinline__
#else
#define VDW_CUT_HD inline
#endif
namespace vdw_cut {
using vdw_next::u32;
struct Counts { unsigned long long words=0,queries=0,rounds=0,input=0,output=0; };
VDW_CUT_HD u32 pop(u32 x){
#if defined(__CUDA_ARCH__)
 return __popc(x);
#else
 return __builtin_popcount(x);
#endif
}
// The two budgets U,V must satisfy U+V <= q-t+1. Every certified left
// mismatch at i<=U and right mismatch at j<=V implies length <=t+i+j-2<q.
// mode=1: balanced budgets; mode=2: left 1, right q-t; mode=3: reverse.
// q<=t gives no pruning. Domain arithmetic must fit signed 32-bit indices;
// the parent geometry and table allocation certify that precondition.
template<class V>
VDW_CUT_HD u32 screen(V v,u32 step,u32 t,u32 q,u32 word,u32 mask,
                      unsigned mode=1,Counts* st=nullptr){
 if(!mask||q<=t||mode==0)return mask;
 if(st){++st->words;st->input+=pop(mask);}
 
 // mode 4 evaluates ALL q-windows containing the anchor. Up to 32 missing
 // positions are handled here; beyond that the mask is conservatively kept.
 // This is exact relative to the available observations, not outside colours.
 if(mode==4){
  const u32 h=q-t;if(h>32){if(st)st->output+=pop(mask);return mask;}
  u32 left[33];left[0]=mask;u32 base_planes[4]{};
  for(u32 k=0;k<v.np;++k)base_planes[k]=v.aligned(k,word);
  auto mismatch=[&](int pos,u32 active)->u32 {
   u32 known=v.valid(pos)&active;if(!known)return 0;
   u32 diff=0;for(u32 k=0;k<v.np;++k)diff|=base_planes[k]^v.window(k,pos);
   if(st)++st->queries;return diff&known;
  };
  int base=(int)(word*32);
  for(u32 j=1;j<=h;++j){left[j]=left[j-1]&~mismatch(base+(int)(j*step),left[j-1]);if(st)++st->rounds;}
  u32 right=mask,out=left[h];
  for(u32 j=1;j<=h&&out!=mask;++j){right&=~mismatch(base-(int)((t-1+j)*step),right);out|=right&left[h-j];if(st)++st->rounds;}
  if(st)st->output+=pop(out);return out;
 }
 const u32 sum=q-t+1;
 u32 U=sum/2, W=sum-U;
 if(mode==2){U=1;W=sum-1;}else if(mode==3){U=sum-1;W=1;}
 const int base=(int)(word*32);
 u32 left=mask,right=mask; // no observed mismatch yet on this side
 u32 planes[4]{};
 for(u32 k=0;k<v.np;++k)planes[k]=v.aligned(k,word);
 const u32 maximum=U>W?U:W;
 for(u32 j=1;j<=maximum&&(left|right);++j){
  if(st)++st->rounds;
  if(j<=U&&left){
   int pos=base+(int)(j*step);u32 known=v.valid(pos)&left;
   if(known){u32 diff=0;for(u32 k=0;k<v.np;++k)diff|=planes[k]^v.window(k,pos);
    left&=~(known&diff);if(st)++st->queries;}
  }
  if(j<=W&&right){
   int pos=base-(int)((t-1+j)*step);u32 known=v.valid(pos)&right;
   if(known){u32 diff=0;for(u32 k=0;k<v.np;++k)diff|=planes[k]^v.window(k,pos);
    right&=~(known&diff);if(st)++st->queries;}
  }
 }
 u32 out=left|right;
 if(st)st->output+=pop(out);
 return out;
}
} // namespace vdw_cut
#undef VDW_CUT_HD
