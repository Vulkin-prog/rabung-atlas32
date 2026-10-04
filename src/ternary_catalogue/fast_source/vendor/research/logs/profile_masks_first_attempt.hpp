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
 VDW_TARGET_HD u32 plane(u32 i)const {u32 s=indices;while(i--)s&=s-1;return __builtin_ctz(s);}
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
template<class C,class F> VDW_TARGET_HD Masks split(C c,F f,u32 selected,u32 step,u32 t,u32 qc,u32 qf,u32 word,u32 anchor,Counters* st=nullptr){
 if(!anchor)return {0,0};
 if(st){++st->words;st->anchor_bits+=vdw_cut::pop(anchor);}
 vdw_cut::Counts cs,fs;
 // Saturation is conservative: it weakens rejection, not coverage.
 if(qc>256)qc=256;if(qf>256)qf=256;
 u32 coarse=vdw_cut::screen(c,step,t,qc,word,anchor,1,st?&cs:nullptr);
 u32 rest=anchor&~coarse;
 Projection<F> extra{f,selected,vdw_cut::pop(selected)};
 u32 fine=refine_anchor(extra,step,t,word,rest,st);
 fine=vdw_cut::screen(f,step,t,qf,word,fine,1,st?&fs:nullptr);
 if(st){st->coarse_bits+=vdw_cut::pop(coarse);st->fine_bits+=vdw_cut::pop(fine);st->coarse_screen_queries+=cs.queries;st->fine_screen_queries+=fs.queries;}
 return {coarse,fine};
}
} // namespace vdw_target
#undef VDW_TARGET_HD
