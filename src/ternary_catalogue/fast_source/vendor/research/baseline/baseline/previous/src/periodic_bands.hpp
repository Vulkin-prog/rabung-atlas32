#pragma once
// Exact rank/select for the PREVIOUS anchor set. No colour-dependent decisions.
#include <cstdint>
#if defined(__CUDACC__)
#define VDW_PHD __host__ __device__ __forceinline__
#else
#define VDW_PHD inline
#endif
namespace vdw_periodic {
using u32=std::uint32_t;using u64=std::uint64_t;using i64=std::int64_t;
struct Range {u32 first=1,last=0;};
struct Task {u32 first=1,last=0,word=0;bool valid=false;};
VDW_PHD u32 words(Range r){return r.first>r.last?0:1+(r.last>>5)-(r.first>>5);}
VDW_PHD i64 floor_div(i64 a,i64 b){i64 q=a/b;return q-((a%b)<0);}
VDW_PHD i64 ceil_div(i64 a,i64 b){return -floor_div(-a,b);}
VDW_PHD u32 gcd(u32 a,u32 b){while(b){u32 c=a%b;a=b;b=c;}return a;}
VDW_PHD u32 period(u32 step){return 32/gcd(step,32);}
VDW_PHD Range clipped(u32 B,u32 minimum,u32 d,u32 g,u32 c,u32 n,u32 A,u32 Z){
 u64 lo=(u64)n*d+A,hi=(u64)n*d+Z;
 if(lo<minimum)lo=minimum;
 if(hi>B)hi=B;
 if(lo>hi||hi<c)return {};
 u32 f=lo>c?(u32)((lo-c+g-1)/g):0,l=(u32)((hi-c)/g);
 return {f,l};
}
// At most two clipped boundary bands; all middle bands are exact translations.
struct Shape {Range begin{},end{};u32 F=0,L=0,step=0,middle=0,bands=0;};
VDW_PHD Shape make_side(u32 B,u32 d,u32 g,u32 c,u32 t,u32 s,u32 phi,u32 A,u32 Z){
 Shape o{};if(!s||!d||!g||A>Z||B<A)return o;
 const i64 flo=ceil_div((i64)A-c,g),fhi=floor_div((i64)Z-c,g);
 if(flo>fhi)return o; // The row has no point even in an unclipped band.
 const u32 minimum=(u32)(((u64)(t-1)*d+1)/2);
 if(minimum>B)return o;
 u32 nlo=minimum>Z?(minimum-Z+d-1)/d:0,nhi=(B-A)/d;
 u32 n0=nlo+(phi+s-nlo%s)%s;if(n0>nhi)return o;
 u32 count=(nhi-n0)/s+1;
 o.begin=clipped(B,minimum,d,g,c,n0,A,Z);o.bands=count;
 if(count==1)return o;
 o.end=clipped(B,minimum,d,g,c,n0+(count-1)*s,A,Z);
 o.middle=count-2;o.F=(u32)((u64)(n0+s)*(d/g)+flo);o.L=(u32)((u64)(n0+s)*(d/g)+fhi);o.step=s*(d/g);
 // No-middle shapes can have hypothetical endpoints outside B. Never accessed.
 return o;
}
struct ShapePair {Shape a{},b{};};
VDW_PHD ShapePair make_shapes(u32 B,u32 d,u32 g,u32 c,u32 t,u32 s,u32 phase){
 ShapePair r{};const u32 minimum=(u32)(((u64)(t-1)*d+1)/2);
 if(minimum>B)return r;
 if(d<=2){if(B>=c){r.a.begin={minimum>c?(minimum-c+g-1)/g:0,(B-c)/g};r.a.bands=1;}return r;}
 const u32 other=(t-2+s-phase)%s;
 if(other==phase)r.a=make_side(B,d,g,c,t,s,phase,1,d-1);
 else{r.a=make_side(B,d,g,c,t,s,phase,1,(d-1)/2);r.b=make_side(B,d,g,c,t,s,other,d/2+1,d-1);}
 return r;
}
// CPU lookup. GPU stores the same <=32 cumulative counts across warp registers.
struct Prepared {Shape shape{};u32 P=1,prefix[33]{},cycle=1,mid_words=0,total=0;};
VDW_PHD u32 band_words(Shape s,u32 k){u64 f=(u64)s.F+(u64)k*s.step,l=(u64)s.L+(u64)k*s.step;return (u32)((l>>5)-(f>>5)+1);}
inline Prepared prepare(Shape s){
 Prepared p{};p.shape=s;
 if(s.middle){p.P=period(s.step);for(u32 k=0;k<p.P;++k)p.prefix[k+1]=p.prefix[k]+band_words(s,k);p.cycle=p.prefix[p.P];p.mid_words=(s.middle/p.P)*p.cycle+p.prefix[s.middle%p.P];}
 else p.prefix[1]=1;
 p.total=words(s.begin)+p.mid_words+words(s.end);return p;
}
inline Task select(const Prepared&p,u32 rank){
 if(rank>=p.total)return {};
 auto s=p.shape;u32 nb=words(s.begin);
 if(rank<nb)return {s.begin.first,s.begin.last,(s.begin.first>>5)+rank,true};
 rank-=nb;
 if(rank<p.mid_words){u32 cycle=rank/p.cycle,rem=rank%p.cycle,k=0;
  for(u32 b=p.P>>1;b;b>>=1)if(p.prefix[k+b]<=rem)k+=b;
  u32 band=cycle*p.P+k,offset=rem-p.prefix[k];
  u32 f=(u32)((u64)s.F+(u64)band*s.step),l=(u32)((u64)s.L+(u64)band*s.step);
  return {f,l,(f>>5)+offset,true};
 }
 rank-=p.mid_words;return {s.end.first,s.end.last,(s.end.first>>5)+rank,true};
}
struct PreparedPair {Prepared a{},b{};u32 total=0;};
inline PreparedPair prepare(ShapePair s){PreparedPair p{prepare(s.a),prepare(s.b),0};p.total=p.a.total+p.b.total;return p;}
inline Task select(const PreparedPair&p,u32 rank){return rank<p.a.total?select(p.a,rank):select(p.b,rank-p.a.total);}
} // namespace
#undef VDW_PHD
