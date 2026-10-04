#pragma once
#include "wide_core.hpp"
#include "../vendor/research/baseline/baseline/previous/src/clipped_cover.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using u32=std::uint32_t;using u64=std::uint64_t;
static inline void require(bool condition,const std::string&s){if(!condition)throw std::runtime_error(s);}
static inline u64 parse64(const char*s){
 require(s&&*s,"empty integer");u64 v=0;for(;*s;++s){require(*s>='0'&&*s<='9',"non-decimal integer");u32 d=*s-'0';require(v<=(UINT64_MAX-d)/10,"integer overflow");v=10*v+d;}return v;
}
static inline u32 parse(const char*s){u64 n=parse64(s);require(n<=UINT32_MAX,"32-bit option overflow");return (u32)n;}
static inline u64 root(u64 p){
 require(p>=1000&&p<(1ull<<34),"wide33 root domain (caller has validated primality)");
 u64 n=p-1;std::vector<u64>factors;for(u64 q=2;q*q<=n;++q)if(n%q==0){factors.push_back(q);do n/=q;while(n%q==0);}if(n>1)factors.push_back(n);
 for(u64 g=2;g<p;++g){bool ok=true;for(u64 q:factors)if(wide33::power128(g,(p-1)/q,p)==1){ok=false;break;}if(ok)return g;}
 throw std::runtime_error("primitive root unavailable");
}
static u32 encoded(u32 colour,u32 r){u32 two=1;while(r%(2*two)==0)two*=2;return colour%two+two*(colour%(r/two));}
static u32 minimal_planes(u32 c,u32 f){
 require(c>=2&&f>c&&f<=9&&f%c==0,"refinement");u32 np=0;while((1u<<np)<f)++np;u32 best=(1u<<np)-1;
 for(u32 s=1;s<(1u<<np);++s){bool ok=true;for(u32 a=0;a<f;++a)for(u32 b=a+1;b<f;++b)if(a%c==b%c&&((encoded(a,f)^encoded(b,f))&s)==0)ok=false;
 if(ok&&(__builtin_popcount(s)<__builtin_popcount(best)||(__builtin_popcount(s)==__builtin_popcount(best)&&s<best)))best=s;}
 return best;
}
static inline void certify_witness33(u64 p,u32 r,u32 maximum,u32 prefix,u64 witness){
 require(p>=1000&&p<(1ull<<34)&&r>=2&&r<=9&&(p-1)%r==0,"witness domain");
 require(maximum&&maximum<=(p-1)/r&&witness&&witness<=(p-1)/2&&((u64)witness+maximum<=p),"witness interval");
 u64 e=(p-1)/r,y=wide33::power128(witness,e,p);
 for(u64 i=1;i<maximum;++i)require(wide33::power128((u64)witness+i,e,p)==y,"witness colour mismatch");
 if(witness>1)require(wide33::power128(witness-1,e,p)!=y,"witness left endpoint");
 if((u64)witness+maximum<p)require(wide33::power128((u64)witness+maximum,e,p)!=y,"witness right endpoint");
 u64 precheck=1;while(precheck+1<p&&wide33::power128(precheck+1,e,p)==1)++precheck;
 require(precheck==prefix,"prefix direct verification");
}
namespace kit {
inline std::vector<std::pair<u32,u32>> groups(u64 p,const std::string &scope){
 if(scope=="24"){require((p-1)%4==0,"24 needs p=1 mod4");return {{2,4}};}
 if(scope=="36"){require((p-1)%6==0,"36 needs p=1 mod6");return {{3,6}};}
 require(scope=="core"||scope=="all","scope must be core, all, 24 or 36");std::vector<std::pair<u32,u32>>v;
 if((p-1)%4==0)v.push_back({2,4});else v.push_back({2,2});if((p-1)%3==0)v.push_back({3,6});
 if(scope=="all")for(u32 r:{5u,7u,8u,9u})if((p-1)%r==0)v.push_back({r,r});return v;
}
inline u32 requested_mask(u64 p,const std::string&s){u32 m=0;for(auto [c,f]:groups(p,s)){m|=1u<<(c-2);m|=1u<<(f-2);}return m;}
inline void dense_host(u64 p,u32 mask,u32*m,u32*pre,u32*witness){
 require(p<=2000000&&wide33::prime(p),"dense CPU oracle domain");u64 g=root(p),x=1;std::vector<unsigned short>idx(p);
 for(u64 i=0;i<p-1;++i){idx[x]=i%2520;x=(unsigned __int128)x*g%p;}
 for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2))){u32 run=0,previous=UINT32_MAX;m[r]=0;pre[r]=1;
 for(u32 j=1;j<p;++j){u32 col=idx[j]%r;run=col==previous?run+1:1;previous=col;if(run>m[r]){m[r]=run;witness[r]=j-run+1;}}
 while((u64)pre[r]+1<p&&idx[pre[r]+1]%r==0)++pre[r];}
}
}

#ifndef KIT_CHARACTER_ROOT
#define KIT_CHARACTER_ROOT 0
#endif
static u64 character_root(u64 p,u32 m){
 require(m>=2&&m<=2520&&(p-1)%m==0,"character root domain");if(m==2)return p-1;
 u32 n=m;std::vector<u32>factors;for(u32 q=2;q*q<=n;++q)if(n%q==0){factors.push_back(q);do n/=q;while(n%q==0);}if(n>1)factors.push_back(n);
 for(u64 a=2;a<p;++a){u64 z=wide33::power128(a,(p-1)/m,p);bool ok=z!=0;for(u32 q:factors)if(wide33::power128(z,m/q,p)==1){ok=false;break;}if(ok)return z;}
 throw std::runtime_error("primitive character root not found");
}
