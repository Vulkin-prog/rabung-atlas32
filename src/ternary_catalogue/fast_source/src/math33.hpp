#if defined(KIT_MONT33_MODE) && KIT_MONT33_MODE !=0
#error wide34 requires the general radix64 Montgomery path
#endif
#pragma once
#ifndef KIT_MONT33_MODE
#define KIT_MONT33_MODE 0
#endif
#include <cstdint>
#if defined(__CUDACC__)
#include <cuda_runtime.h>
#define W33_HD __host__ __device__ __forceinline__
#else
#define W33_HD inline
#endif
namespace wide33 {
using u32=std::uint32_t;using u64=std::uint64_t;using i64=std::int64_t;
inline u64 ninverse(u64 p){u64 x=1;for(int i=0;i<6;++i)x*=2ull-p*x;return 0ull-x;}
constexpr u32 radix_bits=KIT_MONT33_MODE==1?32:64;
// Exact product limbs for canonical 33-bit operands: high limbs are 0 or 1.
W33_HD void product33(u64 a,u64 b,u64 &lo,u64 &hi){
 u64 base=(u64)(u32)a*(u32)b;
 u64 cross=((a>>32)?(u32)b:0ull)+((b>>32)?(u32)a:0ull);
 lo=base+(cross<<32);hi=((a>>32)&(b>>32))+(cross>>32)+(lo<base);
}
// m is arbitrary 64-bit; p is 33-bit. The cross sum cannot overflow u64.
W33_HD u64 product_modulus_hi(u64 m,u64 p){
 u64 base=(u64)(u32)m*(u32)p;
 u64 cross=(u64)(u32)(m>>32)*(u32)p+((p>>32)?(u32)m:0ull);
 u64 lo=base+(cross<<32);
 return ((p>>32)?(m>>32):0ull)+(cross>>32)+(lo<base);
}
// Domain: odd p < 2^33, 0<=a,b<p. The quotient sum is <2p<2^34.
W33_HD u64 mont(u64 a,u64 b,u64 p,u64 ni){
#if KIT_MONT33_MODE==1
 // R=2^32 may be smaller than p. REDC is still exact, with q<3p:
 // (a*b+m*p)/R < p*p/R+p <3p because p<2R. Two subtractions suffice.
 u64 lo,hi;product33(a,b,lo,hi);u32 m=(u32)lo*(u32)ni;
 u64 base=(u64)m*(u32)p,mp=base+((p>>32)?((u64)m<<32):0ull);
 u64 mp_hi=mp<base,sum=lo+mp;
 u64 q=(sum>>32)+((hi+mp_hi+(sum<lo))<<32);
 if(q>=p)q-=p;if(q>=p)q-=p;return q;
#else
 u64 lo=a*b,hi;
#if KIT_MONT33_MODE==2
 product33(a,b,lo,hi);
#elif defined(__CUDA_ARCH__)
 hi=__umul64hi((unsigned long long)a,(unsigned long long)b);
#else
 hi=(unsigned __int128)a*b>>64;
#endif
 u64 m=lo*ni,mhi;
#if KIT_MONT33_MODE==2
 mhi=product_modulus_hi(m,p);
#elif defined(__CUDA_ARCH__)
 mhi=__umul64hi((unsigned long long)m,(unsigned long long)p);
#else
 mhi=(unsigned __int128)m*p>>64;
#endif
 u64 q=hi+mhi+(lo!=0);return q>=p?q-p:q;
#endif
}
// Independent multiplication path. Quotient error <1 for p<2^33;
// the true signed remainder after the estimate lies in (-p,2p).
W33_HD u64 float_mod(u64 a,u64 b,u64 p,double invp){
 u64 q=(u64)((double)a*(double)b*invp);
 i64 rem=(i64)(a*b-q*p);
 if(rem<0)rem+=(i64)p;
 if((u64)rem>=p)rem-=(i64)p;
 return (u64)rem;
}
inline u64 power128(u64 a,u64 e,u64 p){u64 y=1;while(e){if(e&1)y=(unsigned __int128)y*a%p;a=(unsigned __int128)a*a%p;e>>=1;}return y;}
inline bool prime(u64 p){
 if(p<2||p>=(1ull<<34))return false;
 constexpr u64 bases[]={2,3,5,7,11,13,17,19,23,29,31,37};
 for(u64 q:bases)if(p%q==0)return p==q;
 u64 d=p-1;u32 s=0;while(!(d&1)){d>>=1;++s;}
 for(u64 a:bases){u64 x=power128(a,d,p);if(x==1||x==p-1)continue;bool witness=true;
  for(u32 j=1;j<s;++j){x=(unsigned __int128)x*x%p;if(x==p-1){witness=false;break;}}
  if(witness)return false;
 }return true;
}
} // namespace wide33
#undef W33_HD

