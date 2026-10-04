#pragma once
#include "positive_filter.hpp"
#include "sparse_words.hpp"
#if defined(__CUDACC__)
#define VDW_CHD __host__ __device__ __forceinline__
#else
#define VDW_CHD inline
#endif
#ifndef VDW_COMPLETION_MONT
#define VDW_COMPLETION_MONT 0
#endif
#ifndef VDW_COLOUR_FAST
#define VDW_COLOUR_FAST 0
#endif
namespace vdw_next {
using u64=std::uint64_t;using i64=std::int64_t;
struct CompletionView {
 u32 p,r,B,t,D;const std::uint8_t* colour;const u32* inverse;
 u32 ni=0,sigma_value=0xffffffffu;
 VDW_CHD u32 sigma()const{return sigma_value!=0xffffffffu?sigma_value:((p-1)/2)%r;}
 // All colour reductions here have arguments in [0,2r).
 VDW_CHD u32 colour_mod(u32 x)const{
#if VDW_COLOUR_FAST
  return x>=r?x-r:x;
#else
  return x%r;
#endif
 }
 VDW_CHD u32 negative_quotient(u32 b,u32 d)const{
#if VDW_COMPLETION_MONT
  if(ni){
   u64 t=(u64)b*inverse[d],mn=(u64)((u32)t*ni)*p,z=t+mn;
   u32 q=z>>32,y=(z<t||q>=p)?q-p:q;
   return y?p-y:0;
  }
#endif
  u32 y=(u64)b*inverse[d]%p;return y?p-y:0;
 }
 VDW_CHD int small(i64 z)const{
  if(!z||z<-(i64)B||z>B)return -1;
  return colour_mod(colour[z<0?-z:z]+(z<0?sigma():0));
 }
 VDW_CHD int at(u32 x)const{
  if(!x||x>=p)return -1;
  if(x<=B)return colour[x];if(p-x<=B)return colour_mod(colour[p-x]+sigma());
  u32 a=p,b=x;i64 c=0,d=1;
  while(b>B){u32 q=a/b,rem=a-q*b;i64 next=c-(i64)q*d;a=b;b=rem;c=d;d=next;}
  u64 ad=d<0?(u64)-d:(u64)d;
  if(!b||!ad||ad>B)return -1;
  u32 bc=colour_mod(colour[ad]+(d<0?sigma():0));
  return colour_mod(colour[b]+r-bc);
 }
 VDW_CHD int quotient_at(i64 z,u32 x,u32 d)const{
  int c=small(z);return c<0?at(x):colour_mod(c+r-colour[d]);
 }
};
VDW_CHD u64 pack_best(u32 length,u32 start){return ((u64)length<<32)|(0xffffffffu-start);}
struct Completion {u64 best;u32 primitives,heads;bool ok;};
VDW_CHD Completion finish_one(CompletionView v,u64 best,u32 x,i64 a,u32 d){
 Completion result{best,0,0,true};
 if(!x||(u64)x+v.t>v.p){result.ok=false;return result;}
 int c=v.small(a);if(c<0){result.ok=false;return result;}c=v.colour_mod(c+v.r-v.colour[d]);
 if(x>1){int prev=v.quotient_at(a-(i64)d,x-1,d);if(prev<0){result.ok=false;return result;}if(prev==c)return result;}
 result.heads=1;u32 len=v.t;
 while((u64)x+len<v.p){int z=v.quotient_at(a+(i64)len*d,x+len,d);if(z<0){result.ok=false;return result;}if(z!=c)break;++len;}
 u64 b=pack_best(len,x);if(b>result.best)result.best=b;return result;
}
VDW_CHD Completion complete_word(CompletionView v,SparseWord word,u64 seed){
 Completion out{seed,0,0,true};u32 d=word.tagged_d&0x3fffffffu,tag=word.tagged_d>>30;
 u32 g=tag==0?1:tag==1?2:tag==2?3:6;
 if(!d||d>v.D||d%g){out.ok=false;return out;}
 u32 bits=word.left|word.right;
 while(bits){
#if defined(__CUDA_ARCH__)
  u32 k=__ffs(bits)-1;
#else
  u32 k=__builtin_ctz(bits);
#endif
  bits&=bits-1;u32 b=word.first_b+g*k;
  if(!b||b>v.B){out.ok=false;return out;}
  if(!admissible(b,d,v.t))continue;++out.primitives;
  u32 flags=((word.left>>k)&1)|(((word.right>>k)&1)<<1);
  u32 x=v.negative_quotient(b,d);
  if(!x||(u64)x+v.t>v.p){out.ok=false;return out;}
  u32 xr=v.p-x-(v.t-1);
  if(flags&1){auto q=finish_one(v,out.best,x,-(i64)b,d);out.best=q.best;out.heads+=q.heads;if(!q.ok){out.ok=false;return out;}}
  if((flags&2)&&(!(flags&1)||xr!=x)){auto q=finish_one(v,out.best,xr,(i64)b-(v.t-1)*(i64)d,d);out.best=q.best;out.heads+=q.heads;if(!q.ok){out.ok=false;return out;}}
 }
 return out;
}
} // namespace vdw_next
#undef VDW_CHD
