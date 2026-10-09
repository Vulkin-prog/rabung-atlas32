// Separate CPU recheck for the 2026-10-09 bounds addendum.
// No project headers, CUDA, rational covers, or Montgomery arithmetic.
// Primes satisfy 5 <= p < 2^35, r is 2 or 3, and k <= 256.
// Two-bit storage bounds memory; 3 is an unvisited sentinel. Powers use
// exact unsigned 128-bit products. A separate checked bound ensures every
// multiplicative-walk product fits uint64_t. No mirror folding is used.
// The published v1.0 registry and checker are not modified.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>
using U = uint64_t;
static U power(U a, U e, U p) {
    U v=1; using Wide = __uint128_t;
    while(e) {if(e&1) v=static_cast<U>(Wide(v)*a%p); a=static_cast<U>(Wide(a)*a%p); e>>=1;}
    return v;
}
static bool prime(U p) {
    if(p<2) return false;
    if(p%2==0) return p==2;
    for(U q=3;q*q<=p;q+=2) if(p%q==0) return false;
    return true;
}
static std::vector<U> factors(U n) {
    std::vector<U> a;
    for(U q=2;q*q<=n;++q) if(n%q==0) {
        a.push_back(q); do {n/=q;} while(n%q==0);
    }
    if(n>1) a.push_back(n);
    return a;
}
static U integer(const char* s) {
    std::string v(s); size_t end=0;
    if(v.empty() || v[0]<'0' || v[0]>'9') throw std::runtime_error("integer required");
    U n=std::stoull(v,&end);
    if(end!=v.size()) throw std::runtime_error("integer required");
    return n;
}
int main(int argc,char** argv) {
  try {
    if(argc!=4) throw std::runtime_error("usage: two_bit_profile p r k");
    U p=integer(argv[1]), r=integer(argv[2]), k=integer(argv[3]);
    if(p<5 || p>=(U(1)<<35) || r<2 || r>3 || k<3 || k>256 || k>=p)
        throw std::runtime_error("outside checker domain");
    if(!prime(p) || (p-1)%r) throw std::runtime_error("prime/congruence precondition");
    auto fs=factors(p-1);
    U g=2;
    for(;;++g) {
        bool ok=true;
        for(U q:fs) if(power(g,(p-1)/q,p)==1) {ok=false;break;}
        if(ok) break;
    }
    auto begin=std::chrono::steady_clock::now();
    if(g>UINT64_MAX/(p-1)) throw std::runtime_error("walk product capacity");
    std::vector<uint8_t> c(static_cast<size_t>((p+3)/4),255);
    auto get=[&](U y)->unsigned {return (c[static_cast<size_t>(y/4)]>>(2*(y%4)))&3;};
    auto put=[&](U y,unsigned v) {unsigned shift=2*(y%4); auto &b=c[static_cast<size_t>(y/4)];
        b=static_cast<uint8_t>((b & ~(3U<<shift)) | (v<<shift));};
    U x=1, label=0;
    for(U i=0;i<p-1;++i) {
        put(x,static_cast<unsigned>(label));
        if(++label==r) label=0;
        x=x*g%p;
    }
    if(x!=1 || label!=0) throw std::runtime_error("walk closure");
    U M=0,L=0,w=0,run=0,start=1;
    bool prefix=true;
    std::vector<U> counts(static_cast<size_t>(r),0);
    for(U y=1;y<p;++y) {
        unsigned v=get(y);
        if(v>=r) throw std::runtime_error("unvisited residue");
        ++counts[v];
        if(prefix && v==0) ++L; else prefix=false;
        if(y>1 && v==get(y-1)) ++run;
        else {run=1;start=y;}
        if(run>M) {M=run;w=start;}
    }
    for(U v:counts) if(v!=(p-1)/r) throw std::runtime_error("colour population");
    U sigma=((p-1)/2)%r;
    if(sigma!=0 && (r%2 || sigma!=r/2)) throw std::runtime_error("sign class");
    U K=std::max(M+1,(sigma==0?2*L+2:L+2));
    U zeta=power(g,(p-1)/r,p),seed=p;
    std::vector<U> roots(static_cast<size_t>(r),1);
    for(U j=1;j<r;++j) roots[j]=static_cast<U>(__uint128_t(roots[j-1])*zeta%p);
    for(U j=0;j<2048;++j) {
        seed=seed*6364136223846793005ULL+1442695040888963407ULL;
        U y=1+seed%(p-1);
        if(power(y,(p-1)/r,p)!=roots[get(y)])
            throw std::runtime_error("powering comparison");
        if(get(p-y)!=(get(y)+sigma)%r)
            throw std::runtime_error("reflection comparison");
    }
    bool Bstar=true;
    for(U j=0;j<k;++j) {
        bool same=true,have=false; unsigned first=0;
        for(int64_t t=-static_cast<int64_t>(j);t<=static_cast<int64_t>(k-1-j);++t) {
            if(t==0) continue;
            U y=t<0 ? p-static_cast<U>(-t):static_cast<U>(t);
            unsigned v=get(y);
            if(!have) {first=v;have=true;}
            else if(v!=first) same=false;
        }
        if(same) Bstar=false;
    }
    bool accept=M<k && Bstar;
    if(accept!=(K<=k)) throw std::runtime_error("boundary equivalence");
    double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    std::printf("{\"schema\":1,\"p\":%llu,\"r\":%llu,\"k\":%llu,\"g\":%llu,\"M\":%llu,\"L\":%llu,\"K\":%llu,\"w\":%llu,\"Bstar\":%s,\"accepted\":%s,\"strict_bound\":%llu,\"powering_checks\":2048,\"reflection_checks\":2048,\"full_positions\":%llu,\"seconds\":%.6f}\n",
        (unsigned long long)p,(unsigned long long)r,(unsigned long long)k,(unsigned long long)g,
        (unsigned long long)M,(unsigned long long)L,(unsigned long long)K,(unsigned long long)w,
        Bstar?"true":"false",accept?"true":"false",(unsigned long long)((k-1)*p+1),
        (unsigned long long)(p-1),elapsed);
    return accept?0:1;
  } catch(const std::exception& e) {
    std::fprintf(stderr,"REFUSED: %s\n",e.what()); return 2;
  }
}
