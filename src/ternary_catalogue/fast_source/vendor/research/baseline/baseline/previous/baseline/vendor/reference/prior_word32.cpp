// CPU reference of an SM120-oriented word32 filter. NOT a CUDA benchmark.
// Fixed threshold, reflection, exact zero exclusion, early word extinction.
// Build: g++ -O3 -std=c++17 word32_reference.cpp -o word32_reference
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using u32=uint32_t; using u64=uint64_t; using i64=int64_t;
static void require(bool b,const char* s){if(!b)throw std::runtime_error(s);}
static u64 power(u64 a,u64 e,u64 p){u64 z=1;for(;e;e>>=1,a=a*a%p)if(e&1)z=z*a%p;return z;}
static bool prime(u64 p){if(p<2)return false;for(u64 q=2;q*q<=p;++q)if(p%q==0)return false;return true;}
static u32 root(u32 p){u64 n=p-1;std::vector<u64> f;for(u64 q=2;q*q<=n;++q)if(n%q==0){f.push_back(q);while(n%q==0)n/=q;}if(n>1)f.push_back(n);
 for(u32 g=2;g<p;++g){bool ok=true;for(auto q:f)if(power(g,(p-1)/q,p)==1){ok=false;break;}if(ok)return g;}throw std::runtime_error("root");}
struct Table{u32 p,B,D,L,r,g;std::vector<uint8_t> c;std::vector<std::vector<u32>> planes;u32 n;int sigma;
 Table(u32 pp,u32 rr,u32 ll):p(pp),L(ll),r(rr){
  require(prime(p)&&r>=2&&r<=9&&(p-1)%r==0&&L>=3&&L<p&&L<=256,"domain");
  u64 area=(u64)(L-1)*p;B=(u32)std::sqrt((long double)area);
  while((u64)B*B>area){--B;}
  while((u64)(B+1)*(B+1)<=area){++B;}
  ++B;
  D=2*B/(L-1);require(B<p&&D<p,"cover");g=root(p);sigma=((p-1)/2)%r;c.resize(B+1);
  std::vector<u32> spf(B+1);for(u32 q=2;q<=B;++q)if(!spf[q])for(u32 x=q;x<=B;x+=q)if(!spf[x])spf[x]=q;
  std::vector<u64> roots(r);u64 z=power(g,(p-1)/r,p);roots[0]=1;for(u32 j=1;j<r;++j)roots[j]=roots[j-1]*z%p;
  for(u32 x=2;x<=B;++x){u32 q=spf[x];if(q<x)c[x]=(c[q]+c[x/q])%r;else{auto y=power(x,(p-1)/r,p);auto it=std::find(roots.begin(),roots.end(),y);require(it!=roots.end(),"character");c[x]=it-roots.begin();}}
  n=2*B+1;u32 np=0;while((1u<<np)<r)++np;planes.assign(np,std::vector<u32>((n+31)/32+2));
  for(u32 i=0;i<n;++i)if(i!=B){int cc=colour((i64)i-B);for(u32 j=0;j<np;++j)if((cc>>j)&1)planes[j][i>>5]|=1u<<(i&31);}
 }
 int colour(i64 x)const{if(!x||x<-(i64)B||x>B)return -1;return (c[x<0?-x:x]+(x<0?sigma:0))%r;}
 // Equivalent to a 32-bit funnel shift. The shift=0 case avoids C++ UB.
 u32 window(size_t plane,u32 bit)const{u32 w=bit>>5,s=bit&31;u32 a=planes[plane][w];return s?(a>>s)|(planes[plane][w+1]<<(32-s)):a;}
};
struct Result{u32 value=1,seed=1,q=2;bool exact=false;i64 t=1;u32 d=1;
 u64 words=0,stages=0,hits=0,positive_words=0,full_stages=0;};
static std::vector<u32> shifts(u32 q,bool ends){std::vector<u32> v;if(!ends){for(u32 j=1;j<q;++j)v.push_back(j);}else{u32 lo=1,hi=q-1;while(lo<=hi){v.push_back(hi--);if(lo<=hi)v.push_back(lo++);}}return v;}
static u32 candidates(const Table& T,u32 d,u32 q,u32 w,const std::vector<u32>& order,u64& stages,bool poison=false){
 i64 last=2*(i64)T.B-(q-1)*(i64)d;require(last>=0,"nonempty start range");u32 max_index=last/2,start=w*32;
 u32 count=std::min<u32>(32,max_index-start+1),alive=count==32?~0u:((1u<<count)-1);
 if(!poison&&start<=T.B&&T.B<start+32)alive&=~(1u<<(T.B-start));
 std::array<u32,4> base{};for(size_t j=0;j<T.planes.size();++j)base[j]=T.window(j,start);
 for(auto offset:order){if(!alive)break;++stages;u32 bit=start+offset*d,different=0;
  for(size_t j=0;j<T.planes.size();++j)different|=base[j]^T.window(j,bit);
  alive&=~different;
  // Endpoints already constrain the interval; only arithmetic zero is invalid.
  if(!poison&&bit<=T.B&&T.B<bit+32)alive&=~(1u<<(T.B-bit));
 }
 return alive;
}
static Result scan(const Table& T,bool ends=false,bool poison=false){Result R;u32 run=0;int prev=-1;
 for(u32 x=1;x<=T.B;++x){int c=T.colour(x);run=c==prev?run+1:1;prev=c;if(run>R.value){R.value=run;R.t=x-run+1;}}
 R.seed=R.value;R.q=R.seed+1;if(R.value>=T.L)return R;auto order=shifts(R.q,ends);
 for(u32 d=1;d<=T.D;++d){i64 last=2*(i64)T.B-(R.q-1)*(i64)d;if(last<0)continue;
  u32 nw=((u32)(last/2)+32)/32;
  for(u32 w=0;w<nw;++w){++R.words;R.full_stages+=R.q-1;u32 alive=candidates(T,d,R.q,w,order,R.stages,poison);
   if(alive)++R.positive_words;
   while(alive){u32 b=__builtin_ctz(alive);alive&=alive-1;++R.hits;i64 t=32*(i64)w+b-T.B;
    int col=T.colour(t);i64 lo=t,hi=t+(R.q-1)*(i64)d;u32 len=R.q;
    while(len<T.L&&T.colour(lo-d)==col){lo-=d;++len;}
    while(len<T.L&&T.colour(hi+d)==col){hi+=d;++len;}
    if(len>R.value){R.value=len;R.t=lo;R.d=d;}if(R.value>=T.L)return R;
   }
  }
 }
 R.exact=R.value<T.L;return R;
}
static u32 dense(const Table& T){require(T.p<=2000000,"dense reference domain");std::vector<uint8_t> w(T.p);u64 x=1;
 for(u32 j=0;j<T.p-1;++j){w[x]=j%T.r;x=x*T.g%T.p;}require(x==1,"walk");u32 best=0,len=0;int prev=-1;
 for(u32 i=1;i<T.p;++i){len=w[i]==prev?len+1:1;prev=w[i];best=std::max(best,len);}return best;}
static void check(const Table&T,const Result&R,u32 actual){require(R.value<=actual&&std::min(R.value,T.L)==std::min(actual,T.L),"truncation mismatch");
 require(R.exact==(R.value<T.L),"status");u64 inv=power(R.d,T.p-2,T.p),a=(R.t%(i64)T.p+T.p)%T.p,x=a*inv%T.p;
 require(x>=1&&x+R.value<=T.p,"lift domain");u64 z=power(x,(T.p-1)/T.r,T.p);
 for(u32 j=0;j<R.value;++j)require(power(x+j,(T.p-1)/T.r,T.p)==z,"lift colour");}
static u64 test_candidates(const Table&T){u64 checks=0;
 for(u32 q=2;q<=std::min<u32>(9,T.L);++q)for(u32 d=1;d<=T.D;++d){i64 last=2*(i64)T.B-(q-1)*(i64)d;if(last<0)continue;u32 nw=(last/2+32)/32;
  for(u32 w=0;w<nw;++w){u32 expected=0;for(u32 b=0;b<32&&32*w+b<=(u32)(last/2);++b){i64 t=32*w+b-(i64)T.B;int c=T.colour(t);bool ok=c>=0;
    for(u32 j=1;j<q&&ok;++j){ok=T.colour(t+j*(i64)d)==c;}
    if(ok){expected|=1u<<b;}}
   u64 steps=0;require(candidates(T,d,q,w,shifts(q,false),steps)==expected,"word scalar mismatch");
   require(candidates(T,d,q,w,shifts(q,true),steps)==expected,"word permuted mismatch");checks+=2;
  }
 }return checks;}
static void json(const Table&T,const Result&R,bool oracle,u32 actual){std::cout<<"{\"p\":"<<T.p<<",\"r\":"<<T.r<<",\"L\":"<<T.L<<",\"B\":"<<T.B<<",\"D\":"<<T.D
 <<",\"value\":"<<R.value<<",\"seed\":"<<R.seed<<",\"threshold\":"<<R.q<<",\"exact\":"<<(R.exact?"true":"false")
 <<",\"words32\":"<<R.words<<",\"stages32\":"<<R.stages<<",\"unpruned_stages32\":"<<R.full_stages<<",\"positive_words\":"<<R.positive_words<<",\"candidate_bits\":"<<R.hits
 <<",\"dense_oracle\":"<<(oracle?"true":"false")<<",\"oracle_value\":"<<actual<<",\"GPU_benchmark\":false}";}

struct AnchorResult{u32 p,r,t,M,prefix;u64 unique_count,raw_hits,primitive_hits,words,stages,lower,upper;bool dense_oracle;};
static AnchorResult anchors(u32 p,u32 r,const char* output=nullptr){
 require(prime(p)&&r>=2&&r<=9&&(p-1)%r==0,"anchor domain");
 u64 H=(u64)std::sqrt((long double)p);while(H*H<p)++H;while(H>0&&(H-1)*(H-1)>=p)--H;
 u64 Q=1;u32 t=1;u64 lower=0,upper=p-1;
 for(u32 k=2;k<256&&k<p;++k){
  require(Q<=UINT64_MAX/r,"threshold overflow");Q*=r;
  __uint128_t err=(__uint128_t)(Q-1)*(k-1)*H;
  if(err>=p-k){Q/=r;break;}
  t=k;u64 e=(u64)err;lower=(p-k-e+Q-1)/Q;upper=(p-k+e)/Q;
 }
 require(t>=3,"small p: use the dense fallback, threshold < 3");
 Table T(p,r,t);auto order=shifts(t,false);std::vector<u32> inv(T.D+1),starts;
 inv[1]=1;for(u32 d=2;d<=T.D;++d)inv[d]=p-(u64)(p/d)*inv[p%d]%p;
 starts.reserve((size_t)std::min<u64>(2*upper,32000000));
 u64 raw=0,prim=0,words=0,stages=0;
 for(u32 d=1;d<=T.D;++d){i64 last=2*(i64)T.B-(t-1)*(i64)d;if(last<0)continue;u32 nw=(last/2+32)/32;
  for(u32 w=0;w<nw;++w){++words;u32 alive=candidates(T,d,t,w,order,stages);
   while(alive){u32 bit=__builtin_ctz(alive);alive&=alive-1;++raw;i64 a=32*(i64)w+bit-T.B;
    if(std::gcd((u64)(a<0?-a:a),(u64)d)!=1){continue;}
    ++prim;
    u64 num=(a%(i64)p+p)%p,x=num*inv[d]%p;require(x>=1&&x+t<=p,"anchor lift");
    require(starts.size()+2<=32000000,"explicit output memory budget exceeded");
    starts.push_back(x);starts.push_back(p-x-(t-1));
   }
  }
 }
 std::sort(starts.begin(),starts.end());starts.erase(std::unique(starts.begin(),starts.end()),starts.end());
 require(!starts.empty()&&starts.size()>=lower&&starts.size()<=upper,"Weil count/status");
 u32 run=0,best=0,previous=0;for(u32 x:starts){run=(run&&x==previous+1)?run+1:1;best=std::max(best,run);previous=x;}
 u32 M=t-1+best,prefix=0;
 if(starts.front()==1){u32 n=1;while(n<starts.size()&&starts[n]==n+1)++n;prefix=t-1+n;}
 else{while(prefix<T.B&&T.c[prefix+1]==0)++prefix;require(prefix<t,"prefix");}
 if(output){std::ofstream f(output,std::ios::binary);require((bool)f,"anchor file open");
  for(u32 x:starts){char b[4];for(int j=0;j<4;++j)b[j]=(char)((x>>(8*j))&255u);f.write(b,4);}require((bool)f,"anchor file write");}
 bool oracle=p<=2000000;
 if(oracle){require(dense(T)==M,"anchor maximum/dense");std::vector<uint8_t> word(p);u64 x=1;
  for(u32 j=0;j<p-1;++j){word[x]=j%r;x=x*T.g%p;}
  std::vector<u32> want;for(u32 x0=1;x0+t<=p;++x0){bool ok=true;for(u32 j=1;j<t&&ok;++j)ok=word[x0+j]==word[x0];if(ok)want.push_back(x0);}
  require(want==starts,"complete anchor list/dense");u32 pre=1;while(pre+1<p&&word[pre+1]==word[1])++pre;require(pre==prefix,"anchor prefix/dense");
 }
 std::cout<<"{\"status\":\"PASS\",\"p\":"<<p<<",\"r\":"<<r<<",\"anchor_length\":"<<t<<",\"B\":"<<T.B<<",\"D\":"<<T.D
 <<",\"maximum\":"<<M<<",\"prefix\":"<<prefix<<",\"unique_starts\":"<<starts.size()<<",\"raw_hits\":"<<raw<<",\"primitive_hits\":"<<prim
 <<",\"words32\":"<<words<<",\"stages32\":"<<stages<<",\"Weil_lower\":"<<lower<<",\"Weil_upper\":"<<upper
 <<",\"bitplane_bytes\":"<<T.planes.size()*T.planes.front().size()*4
 <<",\"dense_oracle\":"<<(oracle?"true":"false")<<",\"GPU_benchmark\":false}"<<'\n';
 return {p,r,t,M,prefix,starts.size(),raw,prim,words,stages,lower,upper,oracle};
}
int main(int argc,char**argv){try{
 if((argc==4||argc==5)&&std::string(argv[1])=="--anchors"){
  auto parse=[](const char*s){std::string a=s;size_t pos=0;auto v=std::stoull(a,&pos);require(pos==a.size()&&v<=UINT32_MAX,"u32 input");return(u32)v;};
  anchors(parse(argv[2]),parse(argv[3]),argc==5?argv[4]:nullptr);return 0;
 }
 if(argc==2&&std::string(argv[1])=="--test"){u64 profiles=0,sat=0,checks=0,poison_bad=0;
  for(u32 p=5;p<=1000;++p)if(prime(p))for(u32 L:{3u,std::min(8u,p-1),std::min(17u,p-1)})for(u32 r=2;r<=9;++r)if((p-1)%r==0){Table T(p,r,L);u32 a=dense(T);
   for(bool e:{false,true}){auto R=scan(T,e);check(T,R,a);++profiles;sat+=!R.exact;}
   if(p<=100){checks+=test_candidates(T);auto bad=scan(T,false,true);if(bad.value>a)++poison_bad;}
  }
  require(poison_bad>0,"zero poison must be detected");std::cout<<"{\"status\":\"PASS\",\"profile_comparisons\":"<<profiles<<",\"saturations\":"<<sat<<",\"word_mask_comparisons\":"<<checks<<",\"zero_poison_failures_detected\":"<<poison_bad<<",\"GPU_benchmark\":false}\n";return 0;
 }
 if(argc==5&&std::string(argv[1])=="--profile"){auto parse=[](const char*s){std::string a=s;size_t pos=0;auto v=std::stoull(a,&pos);require(pos==a.size()&&v<=UINT32_MAX,"u32 input");return(u32)v;};
  Table T(parse(argv[2]),parse(argv[3]),parse(argv[4]));auto R=scan(T);bool oracle=T.p<=2000000;u32 a=oracle?dense(T):0;if(oracle)check(T,R,a);json(T,R,oracle,a);std::cout<<'\n';return 0;}
 std::cerr<<"Usage: --test | --anchors prime r | --profile prime r L (3<=L<=256, L<prime<2^32)\n";return 2;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
