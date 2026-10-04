#pragma once
// Same integer filter compiled on the CPU and, separately, by nvcc.
// Preconditions are checked by the host. No GPU execution is claimed here.
#include <cstdint>
#if defined(__CUDACC__)
#define VDW_HD __host__ __device__ __forceinline__
#else
#define VDW_HD inline
#endif
namespace vdw {
using u32 = std::uint32_t;
VDW_HD u32 reverse32(u32 x) {
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
VDW_HD u32 low_mask(u32 n) { return n==32 ? ~0u : (n==0 ? 0u : (1u<<n)-1u); }
struct FoldView {
    const u32* data;
    u32 stride, B, np, sign_xor;
    VDW_HD u32 positive(u32 plane,u32 bit) const {
        const u32* p=data+plane*stride;
        const u32 i=bit>>5, shift=bit&31;
        return shift ? (p[i]>>shift)|(p[i+1]<<(32-shift)) : p[i];
    }
    // Bit j represents the encoded colour at a+j; padding/zero are invalid,
    // independently removed by the filter. No shift by 32 is evaluated.
    VDW_HD u32 window(u32 plane,int a) const {
        if(a>=0) return positive(plane,static_cast<u32>(a));
        const u32 flip=(sign_xor&(1u<<plane)) ? ~0u:0u;
        if(a<=-31) return reverse32(positive(plane,static_cast<u32>(-a-31)))^flip;
        const u32 n=static_cast<u32>(-a); // 1..30: the window straddles zero.
        const u32 word=positive(plane,0), neg=low_mask(n);
        return (((reverse32(word)>>(31-n))&neg) | (word<<n)) ^ (flip&neg);
    }
};
VDW_HD u32 equality_keep(u32 alive,u32 a,u32 b) { return alive&~(a^b); }
// A single word of canonical starts a=-B+32*w; t>=3 and d>=1.
VDW_HD u32 folded_candidates(FoldView v,u32 d,u32 t,u32 w,
                            unsigned long long* stages=nullptr) {
    const int last=2*static_cast<int>(v.B)-static_cast<int>((t-1)*d);
    if(last<0) return 0;
    const u32 start=w*32, maximum=static_cast<u32>(last/2);
    if(start>maximum) return 0;
    u32 count=maximum-start+1; if(count>32) count=32;
    u32 alive=low_mask(count), base[4]{};
    const int a=static_cast<int>(start)-static_cast<int>(v.B);
    for(u32 b=0;b<v.np;++b) base[b]=v.window(b,a);
    // Base zero, then every subsequent arithmetic zero, are removed exactly.
    if(a<=0 && a>=-31) alive&=~(1u<<static_cast<u32>(-a));
    for(u32 j=1;j<t && alive;++j) {
        if(stages) ++*stages;
        const int x=a+static_cast<int>(j*d);
        for(u32 b=0;b<v.np;++b) alive=equality_keep(alive,base[b],v.window(b,x));
        if(x<=0 && x>=-31) alive&=~(1u<<static_cast<u32>(-x));
    }
    return alive;
}
} // namespace vdw
#undef VDW_HD
