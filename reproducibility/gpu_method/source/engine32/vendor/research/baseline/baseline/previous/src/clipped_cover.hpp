#pragma once
#include <cstdint>
#include <stdexcept>
namespace vdw_clipped {
using u32=std::uint32_t;using u64=std::uint64_t;
struct Cover{u32 B,D,Q;};
inline bool certified(u32 p,Cover c){
 if(p<3||c.Q<3||c.Q>256||!c.D||!c.B||c.B>=p||c.D>=p||c.B>(1u<<20)||c.D>(1u<<20))return false;
 const u64 H=c.Q-1;
 if(H*c.D>2ull*c.B)return false;
 return 4ull*c.B*c.D-H*c.D*c.D>4ull*p;
}
// Exact integer inequality. No floating-point rounding can certify a bad cover.
inline Cover from_denominator(u32 p,u32 Q,u32 D){
 if(Q<3||Q>256||!D||D>=p||D>(1u<<20))throw std::invalid_argument("clipped cover domain");
 const u64 H=Q-1;
 if(H*D*D>4ull*p)throw std::invalid_argument("denominator beyond lower branch");
 u64 B=(4ull*p+H*D*D)/(4ull*D)+1;
 if(B>=p||B>(1u<<20))throw std::invalid_argument("clipped cover needs dense fallback");
 Cover c{(u32)B,D,Q};if(!certified(p,c))throw std::runtime_error("cover certificate failure");return c;
}
// Minimum D for the given table, on the increasing-volume branch.
inline Cover from_table(u32 p,u32 Q,u32 B){
 if(Q<3||Q>256||!B||B>=p||B>(1u<<20))throw std::invalid_argument("table domain");
 u32 hi=(u32)(2ull*B/(Q-1));if(hi>=p)hi=p-1;
 if(!hi||!certified(p,{B,hi,Q}))throw std::invalid_argument("table too small");
 u32 lo=1;while(lo<hi){u32 mid=lo+(hi-lo)/2;if(certified(p,{B,mid,Q}))hi=mid;else lo=mid+1;}
 return {B,lo,Q};
}
} // namespace
