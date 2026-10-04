#pragma once
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
// Domain: odd p < 2^33, 0<=a,b<p. The quotient sum is <2p<2^34.
W33_HD u64 mont(u64 a,u64 b,u64 p,u64 ni){
 u64 lo=a*b,hi;
#if defined(__CUDA_ARCH__)
 hi=__umul64hi((unsigned long long)a,(unsigned long long)b);
#else
 hi=(unsigned __int128)a*b>>64;
#endif
 u64 m=lo*ni,mhi;
#if defined(__CUDA_ARCH__)
 mhi=__umul64hi((unsigned long long)m,(unsigned long long)p);
#else
 mhi=(unsigned __int128)m*p>>64;
#endif
 // lo + low(m*p) cancels modulo 2^64, with one carry exactly when lo!=0.
 u64 q=hi+mhi+(lo!=0);return q>=p?q-p:q;
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
inline bool prime(u64 p){if(p<2)return false;if(p%2==0)return p==2;if(p%3==0)return p==3;for(u64 d=5;d*d<=p;d+=6)if(p%d==0||p%(d+2)==0)return false;return true;}
} // namespace wide33
#undef W33_HD

