#pragma once
// One candidate may remain useful to the coarse result, to the fine result,
// or to neither. Only the last case is discarded. Calls have no collectives.
#include "../baseline/src/two_cuts.hpp"
#if defined(__CUDACC__)
#define VDW_TARGET_HD __host__ __device__ __forceinline__
#else
#define VDW_TARGET_HD inline
#endif
namespace vdw_target {
using vdw_next::u32;
struct Counters {
 unsigned long long words=0,anchor_bits=0,coarse_bits=0,fine_bits=0;
 unsigned long long residual_queries=0,coarse_screen_queries=0,fine_screen_queries=0;
};
struct Masks {u32 coarse,fine;};
// Projection of the sign-equivariant fine encoding onto selected binary planes.
// Geometry and validity are those of the source row. There are at most 4 planes.
template<class V> struct Projection {
 V source;u32 indices,np;
 VDW_TARGET_HD u32 plane(u32 i)const {u32 s=indices;while(i--)s&=s-1;
#if defined(__CUDA_ARCH__)
 return __ffs(s)-1;
#else
 return __builtin_ctz(s);
#endif
 }
 VDW_TARGET_HD u32 aligned(u32 i,u32 w)const{return source.aligned(plane(i),w);}
 VDW_TARGET_HD u32 positive(u32 i,int b)const{return source.positive(plane(i),b);}
 VDW_TARGET_HD u32 window(u32 i,int b)const{return source.window(plane(i),b);}
 VDW_TARGET_HD u32 valid(int b)const{return source.valid(b);}
};
// Input already passed the entire coarse t-anchor test. Selected fine planes
// need only be injective inside each coarse colour fibre, not on all colours.
template<class V> VDW_TARGET_HD u32 refine_anchor(V v,u32 step,u32 t,u32 word,u32 bits,Counters* st=nullptr){
 if(!bits)return 0;u32 base[4]{};for(u32 k=0;k<v.np;++k)base[k]=v.aligned(k,word);
 const int b=(int)(word*32);
 for(u32 j=1;j<t&&bits;++j){int pos=b-(int)(j*step);
  for(u32 k=0;k<v.np;++k){bits&=~(base[k]^v.window(k,pos));if(st)++st->residual_queries;}
 }
 return bits;
}
// Two nested thresholds, one progressive pass: the strong test continues
// exactly where the weak left/right budgets stopped, reusing registers.
template<class V> VDW_TARGET_HD Masks two_thresholds(V v,u32 step,u32 t,u32 weak,u32 strong,u32 word,u32 anchor,Counters* st=nullptr){
 if(!anchor)return {0,0};
 u32 base[4]{};for(u32 k=0;k<v.np;++k)base[k]=v.aligned(k,word);
 const int b=(int)(word*32);
 u32 left=anchor,right=anchor;
 auto next=[&](u32 &side,int pos){
  u32 known=side&v.valid(pos);if(!known)return;
  u32 diff=0;for(u32 k=0;k<v.np;++k)diff|=base[k]^v.window(k,pos);
  side&=~(known&diff);if(st)++st->coarse_screen_queries;
 };
 u32 U0=0,V0=0;
 if(weak>t){u32 z=weak-t+1;U0=z/2;V0=z-U0;}
 for(u32 j=1;j<=V0&&(left|right);++j){
  if(j<=U0&&left)next(left,b+(int)(j*step));
  if(right)next(right,b-(int)((t-1+j)*step));
 }
 u32 weak_mask=left|right;
 if(!weak_mask||strong<=weak)return {weak_mask,weak_mask};
 u32 z=strong-t+1,U1=z/2,V1=z-U1;
 for(u32 j=1;(U0+j<=U1||V0+j<=V1)&&(left|right);++j){
  if(U0+j<=U1&&left)next(left,b+(int)((U0+j)*step));
  if(V0+j<=V1&&right)next(right,b-(int)((t-1+V0+j)*step));
 }
 return {left|right,weak_mask};
}
template<class C,class F> VDW_TARGET_HD Masks split(C c,F f,u32 selected,u32 step,u32 t,u32 qc,u32 qf,u32 word,u32 anchor,Counters* st=nullptr){
 if(!anchor)return {0,0};
 if(st){++st->words;st->anchor_bits+=vdw_cut::pop(anchor);}
 vdw_cut::Counts fs;
 if(qc>256)qc=256;if(qf>256)qf=256;
 const u32 weak=qc<qf?qc:qf;
 auto nested=two_thresholds(c,step,t,weak,qc,word,anchor,st);
 u32 coarse=nested.coarse;anchor=nested.fine;
 if(!anchor)return {0,0};
 u32 rest=anchor&~coarse;
 Projection<F> extra{f,selected,vdw_cut::pop(selected)};
 u32 fine=refine_anchor(extra,step,t,word,rest,st);
 fine=vdw_cut::screen(f,step,t,qf,word,fine,1,st?&fs:nullptr);
 if(st){st->coarse_bits+=vdw_cut::pop(coarse);st->fine_bits+=vdw_cut::pop(fine);st->fine_screen_queries+=fs.queries;}
 return {coarse,fine};
}
} // namespace vdw_target
#undef VDW_TARGET_HD
