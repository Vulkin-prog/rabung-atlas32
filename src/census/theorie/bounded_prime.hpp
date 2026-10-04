#pragma once
#include <cstdint>
// Deterministic for every uint64_t: the least strong pseudoprime to the
// first twelve prime bases is 318665857834031151167461 > 2^64.
// Sorenson & Webster, arXiv:1509.00864, Theorem 1.1. GPU domain is still <2^33.
inline bool exact_prime64(std::uint64_t n){
 using U=std::uint64_t;
 constexpr U bases[]={2,3,5,7,11,13,17,19,23,29,31,37};
 if(n<2)return false;
 for(U a:bases){if(n==a)return true;if(n%a==0)return false;}
 U d=n-1;unsigned s=0;while(!(d&1)){d>>=1;++s;}
 auto mul=[&](U x,U y)->U{return (unsigned __int128)x*y%n;};
 auto power=[&](U a,U e)->U{U z=1;while(e){if(e&1)z=mul(z,a);a=mul(a,a);e>>=1;}return z;};
 for(U a:bases){U x=power(a,d);if(x==1||x==n-1)continue;
  bool pass=false;for(unsigned j=1;j<s;++j){x=mul(x,x);if(x==n-1){pass=true;break;}if(x==1)return false;}
  if(!pass)return false;
 }return true;
}
