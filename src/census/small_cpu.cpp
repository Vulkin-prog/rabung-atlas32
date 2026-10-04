#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using u64=uint64_t;using u32=uint32_t;
static void need(bool a,const char*s){if(!a)throw std::runtime_error(s);}
static u64 power(u64 a,u64 n,u64 p){u64 v=1;while(n){if(n&1)v=v*a%p;a=a*a%p;n>>=1;}return v;}
static std::vector<u64> factors(u64 n){std::vector<u64> f;for(u64 q=2;q*q<=n;++q)if(n%q==0){f.push_back(q);do n/=q;while(n%q==0);}if(n>1)f.push_back(n);return f;}
static bool prime(u64 p){if(p<2)return false;for(u64 q=2;q*q<=p;++q)if(p%q==0)return false;return true;}
static void emit(u64 p){
 auto a=std::chrono::steady_clock::now();need(p<=2000000&&prime(p),"small CPU domain");
 u64 mask=0,m=1;for(u64 r=2;r<=9;++r)if((p-1)%r==0){mask|=1ull<<(r-2);m=std::lcm(m,r);}
 std::vector<uint16_t> idx(p,65535);
 if(mask){auto f=factors(p-1);u64 g=2;for(;;++g){bool ok=true;for(auto q:f)if(power(g,(p-1)/q,p)==1){ok=false;break;}if(ok)break;need(g<p,"primitive root");}
  u64 x=1;for(u64 e=0;e<p-1;++e){need(idx[x]==65535,"walk duplicate");idx[x]=e%m;x=x*g%p;}need(x==1,"walk end");
 }
 std::ostringstream rows;bool comma=false;
 for(u64 r=2;r<=9;++r)if(mask&(1ull<<(r-2))){
  u64 mx=0,wi=p,pre=0,start=1,len=0;int last=-1;
  auto finish=[&](){if(!len)return;u64 witness=std::min(start,p-(start+len-1));if(len>mx){mx=len;wi=witness;}else if(len==mx)wi=std::min(wi,witness);};
  for(u64 x=1;x<p;++x){int c=idx[x]%r;if(c!=last){finish();start=x;len=1;last=c;}else ++len;}
  finish();while(pre+1<p&&idx[pre+1]%r==0)++pre;
  u64 K=std::max(mx+1,((p-1)/2)%r==0?2*pre+2:pre+2);
  if(comma)rows<<',';comma=true;rows<<"{\"r\":"<<r<<",\"max\":"<<mx<<",\"prefix\":"<<pre<<",\"K\":"<<K<<",\"witness\":"<<wi<<'}';
 }
 double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-a).count();
 std::cout<<"{\"type\":\"measurement\",\"status\":\"PASS\",\"engine\":\"cpu_dense\",\"p\":"<<p<<",\"mask\":"<<mask<<",\"total_ms\":"<<ms<<",\"profiles\":["<<rows.str()<<"]}\n"<<std::flush;
}
int main(){try{
 std::cout<<"{\"type\":\"ready\",\"domain_bits\":0}\n"<<std::flush;std::string line;u64 expected=0;
 while(std::getline(std::cin,line)){if(line=="QUIT"){std::cout<<"{\"type\":\"bye\"}\n"<<std::flush;return 0;}
  std::istringstream in(line);u64 id,p;need(bool(in>>id)&&id==expected,"CPU batch sequence");std::vector<u64> ps;while(in>>p)ps.push_back(p);need(!ps.empty()&&ps.size()<=128,"CPU batch size");auto a=std::chrono::steady_clock::now();for(auto q:ps)emit(q);
  double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-a).count();std::cout<<"{\"type\":\"batch_done\",\"id\":"<<id<<",\"count\":"<<ps.size()<<",\"wall_ms\":"<<ms<<"}\n"<<std::flush;++expected;
 }return 0;
}catch(const std::exception&e){std::cerr<<"VDW_CPU_FAILURE: "<<e.what()<<'\n';return 2;}}
