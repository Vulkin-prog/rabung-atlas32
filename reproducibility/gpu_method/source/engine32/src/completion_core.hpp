#pragma once
// Integration of retained bilateral completion. No global GPU synchronization.
#include "../vendor/research/baseline/baseline/src/coupled_device.hpp"
#if defined(__CUDACC__)
#define KIT_HD __host__ __device__ __forceinline__
#else
#define KIT_HD inline
#endif
namespace kit {
using u32=std::uint32_t; using u64=std::uint64_t;
struct Candidate {u32 tagged_d,first_b,coarse,fine;};
static_assert(sizeof(Candidate)==16,"candidate wire format");
struct State {
 unsigned long long best_c,best_f,snap_c,snap_f,count;
 u32 failed,qc,qf,q,t,cap,D,lo,hi,coarse,fine,epochs_done,prefix_c,prefix_f;
};
KIT_HD bool finish_fine(vdw_next::CompletionView v,vdw_next::SparseWord w,u64 &best){
 u32 d=w.tagged_d&0x3fffffffu,tag=w.tagged_d>>30;
 u32 g=tag==0?1:tag==1?2:tag==2?3:6;
 if(!v.colour||!v.inverse||!d||d>v.D||d>v.B||d%g||w.right)return false;
 u32 bits=w.left;
 while(bits){
#if defined(__CUDA_ARCH__)
  u32 bit=__ffs(bits)-1;
#else
  u32 bit=__builtin_ctz(bits);
#endif
  bits&=bits-1;u64 wide=(u64)w.first_b+(u64)g*bit;
  if(!wide||wide>v.B)return false;u32 b=(u32)wide;
  if(!vdw_next::admissible(b,d,v.t))continue;
  u32 x=v.negative_quotient(b,d);
  if(!x||(u64)x+v.t>v.p)return false;
  std::int64_t a=-(std::int64_t)b;
  int col=v.quotient_at(a,x,d);if(col<0)return false;
  u32 lo=x,hi=x+v.t-1;
  while(lo>1){int z=v.quotient_at(a+((std::int64_t)lo-1-x)*d,lo-1,d);if(z<0)return false;if(z!=col)break;--lo;}
  while(hi<v.p-1){int z=v.quotient_at(a+((std::int64_t)hi+1-x)*d,hi+1,d);if(z<0)return false;if(z!=col)break;++hi;}
  u32 start=lo<v.p-hi?lo:v.p-hi;
  u64 p=vdw_next::pack_best(hi-lo+1,start);if(p>best)best=p;
 }
 return true;
}
KIT_HD bool complete(vdw_next::CompletionView v,u32 coarse,Candidate word,u64 &bc,u64 &bf){
 if(coarse==v.r){if(!finish_fine(v,{word.tagged_d,word.first_b,word.coarse|word.fine,0},bf))return false;bc=bf;return true;}
 if(word.coarse){auto r=vdw_pair::complete(v,coarse,{word.tagged_d,word.first_b,word.coarse,0},bc,bf);if(!r.ok)return false;bc=r.coarse_best;bf=r.fine_best;}
 if(word.fine&&!finish_fine(v,{word.tagged_d,word.first_b,word.fine,0},bf))return false;
 return (bc>>32)>=(bf>>32);
}
KIT_HD void take_snapshot(State &s,u32 lo,u32 hi){
 s.snap_c=s.best_c;s.snap_f=s.best_f;s.qc=s.best_c>>32;s.qf=s.best_f>>32;
 s.q=s.qc<s.qf?s.qc:s.qf;if(s.q>s.cap)s.q=s.cap;
 s.lo=lo;s.hi=hi;s.count=0;
 if(s.q<s.t||lo<1||hi<lo||hi>s.D+1){s.failed=1;if(s.q<s.t)s.q=s.t;}
}
} // namespace kit
#undef KIT_HD
