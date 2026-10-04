#pragma once
// Research candidate. CPU-tested; CUDA qualification is a separate gate.
#include <cstdint>
#if defined(__CUDACC__)
#define VDW_NEW_HD __host__ __device__ __forceinline__
#else
#define VDW_NEW_HD inline
#endif
#ifndef KIT_FUNNEL_LOAD
#define KIT_FUNNEL_LOAD 0
#endif
namespace vdw_next {
using u32=std::uint32_t;
VDW_NEW_HD u32 low(u32 n){return n==32?~0u:(n?((1u<<n)-1):0u);}
VDW_NEW_HD u32 reverse(u32 x){
#if defined(__CUDA_ARCH__)
 return __brev(x);
#else
 x=((x>>1)&0x55555555u)|((x&0x55555555u)<<1);
 x=((x>>2)&0x33333333u)|((x&0x33333333u)<<2);
 x=((x>>4)&0x0f0f0f0fu)|((x&0x0f0f0f0fu)<<4);
 x=((x>>8)&0x00ff00ffu)|((x&0x00ff00ffu)<<8);
 return (x>>16)|(x<<16);
#endif
}
// Each plane has a leading zero word, positive entries 0..B, then two zero words.
// The negative colours follow the validated sign_xor encoding.
struct View {
 const u32* data; u32 stride,B,np,sign_xor;
 VDW_NEW_HD u32 positive(u32 plane,int bit)const {
  // Supports bit >= -31, relying on the explicit leading guard.
  const u32 off=static_cast<u32>(bit+32),s=off&31,i=off>>5;
  const u32* p=data+plane*stride;
  const u32 a=p[i];
#if defined(__CUDA_ARCH__) && KIT_FUNNEL_LOAD
  return __funnelshift_r(a,p[i+1],s);
#else
  return s?(a>>s)|(p[i+1]<<(32-s)):a;
#endif
 }
 VDW_NEW_HD u32 aligned(u32 plane,u32 word)const{return data[plane*stride+1+word];}
 VDW_NEW_HD u32 window(u32 plane,int bit)const {
  if(bit>=0)return positive(plane,bit);
  const u32 flip=(sign_xor&(1u<<plane))?~0u:0u;
  if(bit<=-31)return reverse(positive(plane,-bit-31))^flip;
  const u32 n=static_cast<u32>(-bit),x=positive(plane,0),mask=low(n);
  #if defined(__CUDA_ARCH__) && KIT_FUNNEL_LOAD
  return __funnelshift_l(reverse(x)<<1,x,n)^(flip&mask);
#else
  return (((reverse(x)>>(31-n))&mask)|(x<<n))^(flip&mask);
#endif
 }
 VDW_NEW_HD u32 valid(int bit)const {
  // Bits that refer to nonzero integers in [-B,B]. No invalid read is needed.
  int lo=bit<-(int)B?-(int)B-bit:0;
  int hi=bit+31>(int)B?(int)B-bit+1:32;
  if(lo<0)lo=0;if(hi>32)hi=32;
  if(lo>=hi||hi<=0||lo>=32)return 0;
  u32 m=low((u32)hi)&~low((u32)lo);
  if(bit<=0&&bit>=-31)m&=~(1u<<(u32)(-bit));
  return m;
 }
};
struct Counters {unsigned long long stages=0,direct=0,signed_reads=0,certified_queries=0;};
VDW_NEW_HD u32 domain(View v,u32 d,u32 t,u32 word){
 const u32 first=((t-1)*d+1)/2,b=word*32;
 if(first>v.B||b>v.B||b+31<first)return 0;
 u32 m=low(v.B-b>=31?32:v.B-b+1);
 if(b<first)m&=~low(first-b);
 return m;
}
// Bit k represents b=32*word+k; the canonical old numerator is a=-b.
// The first floor((t-1)/2) comparisons use only positive entries, no reflection.
// Zero exclusion is deliberately deferred to the final primitive predicate.
VDW_NEW_HD u32 filter(View v,u32 d,u32 t,u32 word,Counters* st=nullptr){
 u32 alive=domain(v,d,t,word);if(!alive)return 0;
 u32 base[4]{};for(u32 plane=0;plane<v.np;++plane)base[plane]=v.aligned(plane,word);
 const u32 direct=(t-1)/2;const int b=(int)(32*word);
 for(u32 j=1;j<t&&alive;++j){
  if(st)++st->stages;
  const int x=b-(int)(j*d);
  for(u32 plane=0;plane<v.np;++plane){
   u32 y;if(j<=direct){y=v.positive(plane,x);if(st)++st->direct;}
   else{y=v.window(plane,x);if(st)++st->signed_reads;}
   alive&=~(base[plane]^y);
  }
 }
 return alive;
}
// Equality with base on entries for which the small table gives an exact colour.
// Outside entries are UNKNOWN; callers must keep them, never reject them.
struct Known {u32 mask,equals;};
VDW_NEW_HD Known known(View v,u32 word,int x,u32 alive,Counters*st=nullptr){
 Known k{v.valid(x)&alive,0};if(!k.mask)return k;
 u32 eq=k.mask;
 for(u32 plane=0;plane<v.np;++plane)
  eq&=~(v.aligned(plane,word)^v.window(plane,x));
 k.equals=eq;if(st)++st->certified_queries;return k;
}
struct Oriented {u32 left,right;};
// Each orientation is kept only if it may be a run HEAD with length >= q.
// left: start -b, increasing step d; right: reflected start b-(t-1)*d.
VDW_NEW_HD Oriented prefilter(View v,u32 d,u32 t,u32 q,u32 word,u32 alive,
                              Counters*st=nullptr){
 Oriented out{alive,alive};if(!alive)return out;
 const int b=(int)(32*word),h=(int)t-1;
 // Known equal predecessor disproves headhood. Unknown remains admissible.
 auto kl=known(v,word,b+(int)d,alive,st);
 auto kr=known(v,word,b-(int)(t*d),alive,st);
 out.left&=~kl.equals;out.right&=~kr.equals;
 for(u32 j=t;j<q&&(out.left||out.right);++j){
  if(out.left){auto k=known(v,word,b-(int)(j*d),out.left,st);out.left&=~(k.mask&~k.equals);}
  if(out.right){auto k=known(v,word,b+(int)((j-h)*d),out.right,st);out.right&=~(k.mask&~k.equals);}
 }
 return out;
}
VDW_NEW_HD u32 gcd(u32 a,u32 b){while(b){u32 z=a%b;a=b;b=z;}return a;}
VDW_NEW_HD bool admissible(u32 b,u32 d,u32 t){
 if(gcd(b,d)!=1)return false;
 // With gcd(b,d)=1, a zero term b-j*d can occur only when d=1.
 return d!=1||b>=t;
}
} // namespace vdw_next
#undef VDW_NEW_HD
