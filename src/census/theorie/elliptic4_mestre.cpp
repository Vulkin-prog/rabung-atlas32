#include "bounded_prime.hpp"
#include <algorithm>
#include <numeric>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>
#include <stdexcept>
using U=std::uint64_t;using I=__int128_t;
U powmod(U a,U e,U p){U z=1;while(e){if(e&1)z=(unsigned __int128)z*a%p;a=(unsigned __int128)a*a%p;e>>=1;}return z;}
U isqrt(U n){U x=std::sqrt((long double)n);while((x+1)<=n/(x+1))++x;while(x&&x>n/x)--x;return x;}
struct Point{U x=0,y=0;bool infinity=true;bool operator==(const Point&b)const{return infinity==b.infinity&&(infinity||(x==b.x&&y==b.y));}};
struct Hash{size_t operator()(const Point&a)const{return a.infinity?0:a.x*11400714819323198485ull^(a.y+0x9e3779b97f4a7c15ull);}};
struct Curve{
 U p,s=1;U norm(I x)const{x%=(I)p;if(x<0)x+=p;return x;}
 U inverse(U a)const{I x=0,y=1;U r=p,s=a;while(s){U q=r/s,n=r%s;r=s;s=n;I t=x-(I)q*y;x=y;y=t;}if(r!=1)throw std::runtime_error("inverse");return norm(x);}
 Point negative(Point a)const{if(!a.infinity&&a.y)a.y=p-a.y;return a;}
 Point add(Point a,Point b)const{
  if(a.infinity)return b;if(b.infinity)return a;U slope;
  if(a.x==b.x){if(norm((I)a.y+b.y)==0)return {};if(a.y!=b.y)throw std::runtime_error("curve addition");
   slope=norm((3*(I)a.x*a.x-20*(I)s*a.x+9*(I)s*s)*(I)inverse(norm(2*(I)a.y)));
  }else slope=norm(((I)b.y-a.y)*(I)inverse(norm((I)b.x-a.x)));
  U x=norm((I)slope*slope+10*(I)s-a.x-b.x),y=norm((I)slope*((I)a.x-x)-a.y);return {x,y,false};
 }
 Point multiply(Point a,U n)const{Point z;while(n){if(n&1)z=add(z,a);a=add(a,a);n>>=1;}return z;}
 U rhs(U x)const{return norm((I)x*norm((I)x-s)*norm((I)x-9*(I)s));}
 bool sqrt(U a,U&result)const{
  if(!a){result=0;return true;}if(powmod(a,(p-1)/2,p)!=1)return false;
  if(p%4==3){result=powmod(a,(p+1)/4,p);return true;}
  U q=p-1,s=0;while(!(q&1)){q>>=1;++s;}U z=2;while(powmod(z,(p-1)/2,p)!=p-1)++z;
  U c=powmod(z,q,p),x=powmod(a,(q+1)/2,p),t=powmod(a,q,p),m=s;
  while(t!=1){U i=1,v=(unsigned __int128)t*t%p;while(v!=1&&i<m){v=(unsigned __int128)v*v%p;++i;}if(i==m)throw std::runtime_error("Tonelli");
   U b=powmod(c,1ull<<(m-i-1),p);x=(unsigned __int128)x*b%p;c=(unsigned __int128)b*b%p;t=(unsigned __int128)t*c%p;m=i;
  }result=x;return true;
 }
};
std::vector<U> factors(U n){std::vector<U> v;for(U q=2;q<=n/q;++q)if(n%q==0){v.push_back(q);while(n%q==0)n/=q;}if(n>1)v.push_back(n);return v;}
struct PointOrder{Point point;U order;};
struct Card{U size=0,order=0,attempts=0,annihilator=0,lower=0,upper=0;Point point;U twist_parameter=1;std::vector<PointOrder>witnesses;};
Card cardinality(U p){
 U R=isqrt(4*p),L=p+1-R,H=p+1+R,W=H-L,m=isqrt(W+1);if(m*m<W+1)++m;
 std::mt19937_64 random(p^20261002ull);std::uniform_int_distribution<U> choose(0,p);U twist;
 do{twist=choose(random);}while(!twist||powmod(twist,(p-1)/2,p)!=p-1);
 Curve curves[2]={{p,1},{p,twist}};U orders[2]={1,1};std::vector<PointOrder> witnesses[2];
 for(U attempt=1;;++attempt){U index=(attempt-1)%2;Curve e=curves[index];Point P;
  for(;;){U x=choose(random),y=0;bool sign=random()&1;
   if(x==p){if(!sign){P={};break;}continue;}if(!e.sqrt(e.rhs(x),y)||(!y&&sign))continue;
   P={x,sign&&y?p-y:y,false};break;
  }
  if(P.infinity)continue;
  std::unordered_map<Point,U,Hash> babies;babies.reserve(2*m);Point v;
  for(U j=0;j<m;++j){babies.emplace(v,j);v=e.add(v,P);}Point step=e.negative(e.multiply(P,m)),Q=e.negative(e.multiply(P,L));U n=0;
  for(U i=0;i<=W/m;++i){auto at=babies.find(Q);if(at!=babies.end()&&i*m+at->second<=W){n=L+i*m+at->second;break;}Q=e.add(Q,step);}
  if(!n||!e.multiply(P,n).infinity)throw std::runtime_error("Hasse BSGS");U d=n;
  for(U q:factors(n))while(d%q==0&&e.multiply(P,d/q).infinity)d/=q;
  orders[index]=orders[index]/std::gcd(orders[index],d)*d;witnesses[index].push_back({P,d});
  U first=((L+orders[index]-1)/orders[index])*orders[index];
  if(first<=H&&first+orders[index]>H){Card c;c.size=index?2*p+2-first:first;c.order=orders[index];c.attempts=attempt;c.annihilator=n;c.lower=L;c.upper=H;c.point=P;c.twist_parameter=e.s;c.witnesses=std::move(witnesses[index]);return c;}
 }
}
I chi(U x,U p){if(!(x%p))return 0;return powmod(x%p,(p-1)/2,p)==1?1:-1;}
U quadruples(U p,U N){I boundary=4*chi(2,p)+2*chi(6,p)+2*chi(3,p)+2*chi(p-1,p)+2*chi(p-2,p);
 I S=(I)N-p-2,numerator=(I)p-10-boundary+S;
 if(numerator<0||numerator%8)throw std::runtime_error("quadratic moment congruence");return numerator/8;
}
void print(U p,const Card&c,double seconds){
 std::cout<<"{\"status\":\"PASS\",\"p\":"<<p<<",\"r\":2,\"quadruples\":"<<quadruples(p,c.size)<<",\"elliptic_points\":"<<c.size<<",\"order_lcm\":"<<c.order<<",\"hasse_lower\":"<<c.lower<<",\"hasse_upper\":"<<c.upper<<",\"attempts\":"<<c.attempts<<",\"twist_parameter\":"<<c.twist_parameter<<",\"certificate_points\":[";
 bool first=true;for(auto w:c.witnesses){if(!first)std::cout<<',';first=false;std::cout<<"{\"x\":"<<w.point.x<<",\"y\":"<<w.point.y<<",\"order\":"<<w.order<<",\"prime_factors\":[";bool f=true;
  for(U q:factors(w.order)){if(!f)std::cout<<',';f=false;std::cout<<q;}std::cout<<"]}";}
 std::cout<<"],\"certificate\":\"unique_lcm_multiple_in_Hasse\",\"seconds\":"<<seconds<<"}\n";
}
int main(int argc,char**argv){try{
 if(argc>1){for(int i=1;i<argc;++i){U p=std::stoull(argv[i]);if(p<1000||p>=(1ull<<33)||!exact_prime64(p))throw std::runtime_error("point-count domain");auto begin=std::chrono::steady_clock::now();auto c=cardinality(p);print(p,c,std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count());}return 0;}
 auto begin=std::chrono::steady_clock::now();U tested=0,fallbacks=0,entries=0;
 std::vector<unsigned char> sieve(200001,1);sieve[0]=sieve[1]=0;for(U q=2;q*q<=200000;++q)if(sieve[q])for(U x=q*q;x<=200000;x+=q)sieve[x]=0;
 for(U p=1009;p<=200000;++p)if(sieve[p]&&(p<=5000||p%2520==1)){
  auto c=cardinality(p);fallbacks+=!c.order;Curve e{p};I sum=0,S4=0;U A4=0;
  std::vector<I> character(p);for(U x=1;x<p;++x)character[x]=chi(x,p);
  for(U x=0;x<p;++x){sum+=character[e.rhs(x)];U product=1;for(U j=0;j<4;++j)product=(unsigned __int128)product*((x+j)%p)%p;S4+=character[product];}
  for(U x=1;x+3<p;++x)A4+=character[x]==character[x+1]&&character[x]==character[x+2]&&character[x]==character[x+3];
  if(c.size!=(I)p+1+sum||c.size!=(I)p+2+S4||A4!=quadruples(p,c.size))throw std::runtime_error("dense fourth-moment mismatch");
  Curve proof{p,c.twist_parameter};for(auto w:c.witnesses){
   if(proof.rhs(w.point.x)!=(unsigned __int128)w.point.y*w.point.y%p||!proof.multiply(w.point,w.order).infinity)throw std::runtime_error("certificate");
   for(U q:factors(w.order))if(proof.multiply(w.point,w.order/q).infinity)throw std::runtime_error("point order");
  }
  ++tested;entries+=p;
 }
 std::cout<<"{\"status\":\"PASS\",\"primes\":"<<tested<<",\"dense_entries\":"<<entries<<",\"small_dense_fallbacks\":"<<fallbacks<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<"}\n";
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
