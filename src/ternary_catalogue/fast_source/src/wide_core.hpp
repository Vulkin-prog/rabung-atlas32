#pragma once
#include "math33.hpp"
#include "../vendor/research/src/profile_masks.hpp"
#include "../vendor/research/baseline/baseline/previous/baseline/vendor/src/sparse_words.hpp"
#if defined(__CUDACC__)
#define KIT_HD __host__ __device__ __forceinline__
#else
#define KIT_HD inline
#endif
namespace vdw_next {
using u64=std::uint64_t;using i64=std::int64_t;
struct CompletionView {
 u64 p;u32 r,B,t,D;const std::uint8_t*colour;const u64*inverse;
 u64 ni=0;u32 sigma_value=0xffffffffu;
 KIT_HD u32 sigma()const{return sigma_value!=0xffffffffu?sigma_value:(u32)(((p-1)/2)%r);}
 KIT_HD u32 colour_mod(u32 x)const{return x>=r?x-r:x;}
 KIT_HD u64 negative_quotient(u32 b,u32 d)const{u64 y=wide33::mont(b,inverse[d],p,ni);return y?p-y:0;}
 KIT_HD int small(i64 z)const{if(!z||z<-(i64)B||z>B)return -1;return colour_mod(colour[z<0?-z:z]+(z<0?sigma():0));}
 KIT_HD int at(u64 x)const{
  if(!x||x>=p)return -1;if(x<=B)return colour[x];if(p-x<=B)return colour_mod(colour[p-x]+sigma());
  u64 a=p,b=x;i64 c=0,d=1;
  while(b>B){u64 q=a/b,rem=a-q*b;i64 next=c-(i64)q*d;a=b;b=rem;c=d;d=next;}
  u64 ad=d<0?(u64)-d:(u64)d;if(!b||!ad||ad>B)return -1;
  u32 bc=colour_mod(colour[ad]+(d<0?sigma():0));return colour_mod(colour[b]+r-bc);
 }
 KIT_HD int quotient_at(i64 z,u64 x,u32 d)const{int c=small(z);return c<0?at(x):colour_mod(c+r-colour[d]);}
};
// Qualified wide34 representation:16 length bits and48 witness bits.
// Runs beyond65535 or invalid intervals cause refusal, never truncation.
KIT_HD u64 pack_best(u32 length,u64 start){return ((u64)length<<48)|(0xffffffffffffull-start);}
}
namespace kit {
using u32=std::uint32_t;using u64=std::uint64_t;using i64=std::int64_t;
struct Candidate {u32 tagged_d,first_b,coarse,fine;};
static_assert(sizeof(Candidate)==16,"candidate wire format");
struct State {
 unsigned long long best_c,best_f,snap_c,snap_f,count;
 u32 failed,qc,qf,q,t,cap,D,lo,hi,coarse,fine,epochs_done,prefix_c,prefix_f;
 u32 stop_k,stop_r,stop_ab,stop_flag,stop_fusion,fixed_q;
};

KIT_HD bool packed_run(vdw_next::CompletionView v,u64 lo,u64 hi,u64 &best){
 if(!lo||hi<lo||hi>=v.p)return false;u64 len=hi-lo+1,mirror=v.p-hi,start=lo<mirror?lo:mirror;
 if(len>65535ull||len>(v.p-1)/2||!start||start>0xffffffffffffull)return false;
 u64 z=vdw_next::pack_best((u32)len,start);if(z>best)best=z;return true;
}
KIT_HD bool finish_fine(vdw_next::CompletionView v,vdw_next::SparseWord w,u64 &best){
 u32 d=w.tagged_d&0x3fffffffu,tag=w.tagged_d>>30,g=tag==0?1:tag==1?2:tag==2?3:6;
 if(!v.colour||!v.inverse||!d||d>v.D||d>v.B||d%g||w.right)return false;u32 bits=w.left;
 while(bits){
#if defined(__CUDA_ARCH__)
  u32 bit=__ffs(bits)-1;
#else
  u32 bit=__builtin_ctz(bits);
#endif
  bits&=bits-1;u64 wide=(u64)w.first_b+(u64)g*bit;if(!wide||wide>v.B)return false;u32 b=(u32)wide;
  if(!vdw_next::admissible(b,d,v.t))continue;u64 x=v.negative_quotient(b,d);if(!x||x+v.t>v.p)return false;
  i64 a=-(i64)b;int col=v.quotient_at(a,x,d);if(col<0)return false;u64 lo=x,hi=x+v.t-1;
  while(lo>1){int z=v.quotient_at(a+((i64)lo-1-(i64)x)*d,lo-1,d);if(z<0)return false;if(z!=col)break;--lo;}
  while(hi<v.p-1){int z=v.quotient_at(a+((i64)hi+1-(i64)x)*d,hi+1,d);if(z<0)return false;if(z!=col)break;++hi;}
  if(!packed_run(v,lo,hi,best))return false;
 }
 return true;
}
KIT_HD int coarse_colour(int x,u32 c){return x>=(int)c?x-(int)c:x;}
KIT_HD bool complete(vdw_next::CompletionView v,u32 coarse,Candidate word,u64 &bc,u64 &bf){
 if(v.p>=(1ull<<34)||coarse<2||v.r>9||v.r<coarse||v.r%coarse||v.t<3||v.t>=v.p)return false;
 if(coarse==v.r){if(!finish_fine(v,{word.tagged_d,word.first_b,word.coarse|word.fine,0},bf))return false;bc=bf;return true;}
 // Current paired groups are exactly f=2*c; the colour reduction uses that bound.
 if(v.r!=2*coarse)return false;
 u32 d=word.tagged_d&0x3fffffffu,tag=word.tagged_d>>30,g=tag==0?1:tag==1?2:tag==2?3:6;
 if(!v.colour||!v.inverse||!d||d>v.D||d>v.B||d%g)return false;
 u32 bits=word.coarse;
 while(bits){
#if defined(__CUDA_ARCH__)
  u32 k=__ffs(bits)-1;
#else
  u32 k=__builtin_ctz(bits);
#endif
  bits&=bits-1;u64 wide=(u64)word.first_b+(u64)g*k;if(!wide||wide>v.B)return false;u32 b=(u32)wide;
  if(!vdw_next::admissible(b,d,v.t))continue;u64 x=v.negative_quotient(b,d);if(!x||x+v.t>v.p)return false;
  i64 a=-(i64)b;int cf=v.quotient_at(a,x,d);if(cf<0)return false;int colour=coarse_colour(cf,coarse);u64 l=x,h=x+v.t-1;
  while(l>1){int c=v.quotient_at(a+((i64)l-1-(i64)x)*d,l-1,d);if(c<0)return false;if(coarse_colour(c,coarse)!=colour)break;--l;}
  while(h<v.p-1){int c=v.quotient_at(a+((i64)h+1-(i64)x)*d,h+1,d);if(c<0)return false;if(coarse_colour(c,coarse)!=colour)break;++h;}
  if(!packed_run(v,l,h,bc))return false;if(h-l+1<(bf>>48))continue;
  int previous=-1;u32 run=0;
  for(u64 pos=l;pos<=h;++pos){int c=v.quotient_at(a+((i64)pos-(i64)x)*d,pos,d);if(c<0)return false;
   if(c==previous){if(run==0xffffffffu)return false;++run;}else run=1;previous=c;
   if(!packed_run(v,pos-run+1,pos,bf))return false;
  }
 }
 if(word.fine&&!finish_fine(v,{word.tagged_d,word.first_b,word.fine,0},bf))return false;
 return (bc>>48)>=(bf>>48);
}
KIT_HD void take_snapshot(State &s,u32 lo,u32 hi){
 s.snap_c=s.best_c;s.snap_f=s.best_f;s.qc=s.best_c>>48;s.qf=s.best_f>>48;
 s.q=s.qc<s.qf?s.qc:s.qf;if(s.q>s.cap)s.q=s.cap;
 if(s.fixed_q){s.qc=s.qf=s.q=s.fixed_q;} // Filter lengths, never fictitious best-run witnesses.
 s.lo=lo;s.hi=hi;s.count=0;
 if(s.q<s.t||lo<1||hi<lo||hi>s.D+1){s.failed=1;if(s.q<s.t)s.q=s.t;}
}
} // namespace kit
#undef KIT_HD
