#pragma once
// A geometric hitting set of short windows inside every longer window.
// No character values or probabilistic hypothesis are used by this planner.
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>
namespace vdw_anchor {
using u32=std::uint32_t;using u64=std::uint64_t;
struct Band {u32 first,last;};
struct Plan {u32 d,phase,spacing;u64 word_cost;};
inline bool selected(std::int64_t b,u32 d,u32 t,u32 s,u32 phase){
 if(d<=2)return true;
 auto n=b/(std::int64_t)d,c=b%(std::int64_t)d;if(c<0){c+=d;--n;}
 if(!c||2*c==d)return false;
 u32 phi=2*c<d?phase:(t-2+s-(phase%s))%s;
 auto z=n%(std::int64_t)s;if(z<0)z+=s;
 return z==phi;
}
// Canonical b >= (t-1)d/2. Reflected anchors h*d-b have the same selection.
inline std::vector<Band> bands(u32 B,u32 d,u32 t,u32 s,u32 phase){
 if(!d||t<3||!s||phase>=s)throw std::invalid_argument("anchor band parameters");
 const u64 lower=((u64)(t-1)*d+1)/2;
 std::vector<Band> out;if(lower>B)return out;
 if(d<=2){out.push_back({(u32)lower,B});return out;}
 const u32 other=(t-2+s-phase)%s;
 auto add=[&](u64 lo,u64 hi){lo=std::max<u64>(lo,lower);hi=std::min<u64>(hi,B);if(lo<=hi)out.push_back({(u32)lo,(u32)hi});};
 // Increasing n keeps the bands disjoint and ordered. d divides no accepted b.
 for(u32 n=(u32)(lower/d);n<=B/d;++n){
  u64 a=(u64)n*d;
  if(phase==other){if(n%s==phase)add(a+1,a+d-1);}
  else{
   if(n%s==phase)add(a+1,a+(d-1)/2);
   if(n%s==other)add(a+d/2+1,a+d-1);
  }
 }
 return out;
}
inline u64 word_cost(const std::vector<Band>& v,u32 d,bool wheel){
 const u32 g=wheel?std::gcd(d,6u):1;u64 cost=0;
 for(auto band:v)for(u32 c=(g==1?0:1);c<g;++c){
  if(g>1&&std::gcd(c,g)!=1)continue;
  if(band.last<c)continue;
  u32 a=band.first>c?(band.first-c+g-1)/g:0,b=(band.last-c)/g;
  if(a<=b)cost+=(b>>5)-(a>>5)+1;
 }
 return cost;
}
inline Plan choose_reference(u32 B,u32 d,u32 t,u32 s,bool wheel,bool optimize){
 Plan p{d,0,s,~u64(0)};
 for(u32 phase=0;phase<(optimize?s:1);++phase){
  u64 c=word_cost(bands(B,d,t,s,phase),d,wheel);
  if(c<p.word_cost){p.phase=phase;p.word_cost=c;}
 }
 return p;
}
inline u64 one_band_cost(u64 lo,u64 hi,u32 B,u32 lower,u32 d,bool wheel){
 lo=std::max<u64>(lo,lower);hi=std::min<u64>(hi,B);if(lo>hi)return 0;
 u32 g=wheel?std::gcd(d,6u):1;u64 cost=0;
 for(u32 c=g==1?0:1;c<g;++c){if(g>1&&std::gcd(c,g)!=1)continue;if(hi<c)continue;
  u64 a=lo>c?(lo-c+g-1)/g:0,b=(hi-c)/g;if(a<=b)cost+=(b>>5)-(a>>5)+1;
 }
 return cost;
}
// Accumulate all phases in ONE pass: each half-band contributes to one phase.
inline Plan choose(u32 B,u32 d,u32 t,u32 s,bool wheel,bool optimize){
 if(!optimize||d<=2)return choose_reference(B,d,t,s,wheel,false);
 std::vector<u64> costs(s,0);u32 lower=(u32)(((u64)(t-1)*d+1)/2);
 for(u32 n=lower/d;n<=B/d;++n){u32 a=n%s,b=(t-2+s-a)%s;u64 start=(u64)n*d;
  if(a==b)costs[a]+=one_band_cost(start+1,start+d-1,B,lower,d,wheel);
  else{costs[a]+=one_band_cost(start+1,start+(d-1)/2,B,lower,d,wheel);
       costs[b]+=one_band_cost(start+d/2+1,start+d-1,B,lower,d,wheel);}
 }
 Plan best{d,0,s,costs[0]};for(u32 a=1;a<s;++a)if(costs[a]<best.word_cost){best.phase=a;best.word_cost=costs[a];}
 return best;
}
} // namespace vdw_anchor
