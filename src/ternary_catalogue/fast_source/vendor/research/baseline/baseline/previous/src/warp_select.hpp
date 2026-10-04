#pragma once
#include "periodic_bands.hpp"
#if defined(__CUDACC__)
#define VDW_DHD __host__ __device__ __forceinline__
#else
#define VDW_DHD inline
#endif
#ifndef KIT_RECIP_SELECT
#define KIT_RECIP_SELECT 0
#endif
namespace vdw_periodic {
struct DevicePlan{Shape shape{};u32 P=1,cycle=1,mid_words=0,total=0;float inverse_cycle=0;};
// All threads of a warp call this with the SAME plan and possibly different ranks.
// Fetch(src) shuffles one cumulative count from src; all fetch calls are unconditional.
// No array of 33 prefixes is stored per CUDA thread.
template<class Fetch>
VDW_DHD Task warp_select(DevicePlan p,u32 rank,Fetch fetch){
 const bool valid=rank<p.total;const u32 nb=words(p.shape.begin);
 const bool middle=valid&&rank>=nb&&rank-nb<p.mid_words;
 const u32 mr=middle?rank-nb:0;u32 cycle,rem;
#if defined(__CUDA_ARCH__) && KIT_RECIP_SELECT
 if(p.inverse_cycle){
  cycle=__float2uint_rz(__fmul_rn((float)mr,p.inverse_cycle));
  u32 product=cycle*p.cycle;
  if(product>mr){--cycle;product-=p.cycle;}
  rem=mr-product;if(rem>=p.cycle){++cycle;rem-=p.cycle;}
 }else{cycle=mr/p.cycle;rem=mr%p.cycle;}
#else
 cycle=mr/p.cycle;rem=mr%p.cycle;
#endif
 u32 k=0;
 for(u32 b=p.P>>1;b;b>>=1){u32 end=fetch(k+b-1);if(end<=rem)k+=b;}
 u32 previous=fetch(k?k-1:0);if(!k)previous=0;
 if(!valid)return {};
 if(rank<nb)return {p.shape.begin.first,p.shape.begin.last,(p.shape.begin.first>>5)+rank,true};
 if(middle){u32 band=cycle*p.P+k,f=p.shape.F+band*p.shape.step,l=p.shape.L+band*p.shape.step;return {f,l,(f>>5)+(rem-previous),true};}
 u32 offset=rank-nb-p.mid_words;return {p.shape.end.first,p.shape.end.last,(p.shape.end.first>>5)+offset,true};
}
} // namespace
#undef VDW_DHD
