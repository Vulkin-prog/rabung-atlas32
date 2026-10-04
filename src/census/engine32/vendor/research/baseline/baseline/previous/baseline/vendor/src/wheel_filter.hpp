#pragma once
#ifndef KIT_FILTER_MODE
#define KIT_FILTER_MODE 0
#endif
#include "positive_filter.hpp"
#if defined(__CUDACC__)
#define VDW_WHD __host__ __device__ __forceinline__
#else
#define VDW_WHD inline
#endif
namespace vdw_next {
// For g>1, c is a unit residue 1..g-1. Row n encodes C(g*n+c).
// Negative n uses the positive row g-c. There is no zero integer in this row.
struct ResidueView {
 const u32 *data,*mirror;u32 stride,np,sign_xor,maximum,mirror_maximum;
 VDW_WHD u32 read(const u32*p,u32 plane,int bit)const{
  u32 off=(u32)(bit+32),i=off>>5,s=off&31;const u32*z=p+plane*stride;
  u32 a=z[i];
#if defined(__CUDA_ARCH__) && KIT_FUNNEL_LOAD
  return __funnelshift_r(a,z[i+1],s);
#else
  return s?(a>>s)|(z[i+1]<<(32-s)):a;
#endif
 }
 VDW_WHD u32 positive(u32 plane,int bit)const{return read(data,plane,bit);}
 VDW_WHD u32 aligned(u32 plane,u32 word)const{return data[plane*stride+1+word];}
 VDW_WHD u32 window(u32 plane,int bit)const{
  if(bit>=0)return positive(plane,bit);
  u32 flip=(sign_xor&(1u<<plane))?~0u:0u;
  if(bit<=-32)return reverse(read(mirror,plane,-bit-32))^flip;
  u32 n=(u32)-bit,m=low(n),a=read(mirror,plane,0),b=positive(plane,0);
  #if defined(__CUDA_ARCH__) && KIT_FUNNEL_LOAD
  return __funnelshift_l(reverse(a),b,n)^(flip&m);
#else
  return ((reverse(a)>>(32-n))|(b<<n))^(flip&m);
#endif
 }
 VDW_WHD u32 valid(int bit)const{
  const int minimum=-(int)mirror_maximum-1;
  int lo=bit<minimum?minimum-bit:0;
  int hi=bit+31>(int)maximum?(int)maximum-bit+1:32;
  if(lo<0)lo=0;if(hi>32)hi=32;
  if(lo>=hi||hi<=0||lo>=32)return 0;
  return low((u32)hi)&~low((u32)lo);
 }
};
VDW_WHD u32 range_mask(u32 first,u32 maximum,u32 word){
 u32 x=word*32;if(first>maximum||x>maximum||x+31<first)return 0;
 u32 m=low(maximum-x>=31?32:maximum-x+1);if(x<first)m&=~low(first-x);return m;
}
template<class V>
VDW_WHD u32 filter_range(V v,u32 d,u32 t,u32 word,u32 first,u32 maximum,Counters*st=nullptr){
 u32 alive=range_mask(first,maximum,word);if(!alive)return 0;
 u32 base[4]{};for(u32 p=0;p<v.np;++p)base[p]=v.aligned(p,word);
 int b=(int)(32*word);u32 direct=(t-1)/2;
#if KIT_FILTER_MODE==1 || KIT_FILTER_MODE==2 || KIT_FILTER_MODE==4 || KIT_FILTER_MODE==5
#if KIT_FILTER_MODE>=4
 #pragma unroll 2
#endif
 for(u32 j=1;j<=direct&&alive;++j){if(st)++st->stages;int x=b-(int)(j*d);
  for(u32 p=0;p<v.np
#if KIT_FILTER_MODE==2 || KIT_FILTER_MODE==5
   &&alive
#endif
  ;++p){alive&=~(base[p]^v.positive(p,x));if(st)++st->direct;}
 }
#if KIT_FILTER_MODE>=4
 #pragma unroll 2
#endif
 for(u32 j=direct+1;j<t&&alive;++j){if(st)++st->stages;int x=b-(int)(j*d);
  for(u32 p=0;p<v.np
#if KIT_FILTER_MODE==2 || KIT_FILTER_MODE==5
   &&alive
#endif
  ;++p){alive&=~(base[p]^v.window(p,x));if(st)++st->signed_reads;}
 }
#else
 for(u32 j=1;j<t&&alive;++j){if(st)++st->stages;int x=b-(int)(j*d);
  for(u32 p=0;p<v.np
#if KIT_FILTER_MODE==3
   &&alive
#endif
  ;++p){u32 y;
   if(j<=direct){y=v.positive(p,x);if(st)++st->direct;}
   else{y=v.window(p,x);if(st)++st->signed_reads;}
   alive&=~(base[p]^y);
  }
 }
#endif
 return alive;
}
template<class V>
VDW_WHD Known known_row(V v,u32 word,int x,u32 alive,Counters*st=nullptr){
 Known k{v.valid(x)&alive,0};if(!k.mask)return k;u32 eq=k.mask;
 for(u32 p=0;p<v.np;++p)eq&=~(v.aligned(p,word)^v.window(p,x));
 k.equals=eq;if(st)++st->certified_queries;return k;
}
template<class V>
VDW_WHD Oriented prefilter_row(V v,u32 d,u32 t,u32 q,u32 word,u32 bits,Counters*st=nullptr){
 Oriented out{bits,bits};if(!bits)return out;int b=(int)(32*word),h=(int)t-1;
 auto k=known_row(v,word,b+(int)d,bits,st);out.left&=~k.equals;
 k=known_row(v,word,b-(int)(t*d),bits,st);out.right&=~k.equals;
 for(u32 j=t;j<q&&(out.left||out.right);++j){
  if(out.left){k=known_row(v,word,b-(int)(j*d),out.left,st);out.left&=~(k.mask&~k.equals);}
  if(out.right){k=known_row(v,word,b+(int)((j-h)*d),out.right,st);out.right&=~(k.mask&~k.equals);}
 }
 return out;
}
} // namespace vdw_next
#undef VDW_WHD
