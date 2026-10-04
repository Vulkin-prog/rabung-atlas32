#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>
using i64=int64_t;using u64=uint64_t;
static constexpr int N=1000000,A=100,X=100000;
static std::vector<int> ps,pi;static std::vector<int> phi_table;static std::unordered_map<u64,i64> memo;
static u64 root(u64 x,int k){u64 r=std::pow((long double)x,1.L/k);auto less=[&](u64 a){__uint128_t p=1;for(int j=0;j<k;++j)p*=a;return p<=x;};while(less(r+1))++r;while(!less(r))--r;return r;}
static i64 phi(u64 x,int a){if(!a)return x;if(a<A&&x<X)return phi_table[a*X+x];return phi(x,a-1)-phi(x/ps[a-1],a-1);}
static i64 count(u64 x){
 if(x<N)return pi[x];auto found=memo.find(x);if(found!=memo.end())return found->second;
 int a=count(root(x,4)),b=count(root(x,2)),c=count(root(x,3));
 i64 total=phi(x,a)+(i64)(b+a-2)*(b-a+1)/2;
 for(int i=a;i<b;++i){u64 w=x/ps[i];total-=count(w);if(i<c){int lim=count(root(w,2));for(int j=i;j<lim;++j)total-=count(w/ps[j])-j;}}
 return memo[x]=total;
}
int main(int argc,char**argv){
 std::vector<bool> mark(N,true);mark[0]=mark[1]=false;pi.resize(N);
 for(int i=2;i<N;++i){if(mark[i]){ps.push_back(i);if((i64)i*i<N)for(int j=i*i;j<N;j+=i)mark[j]=false;}pi[i]=ps.size();}
 phi_table.resize(A*X);for(int x=0;x<X;++x)phi_table[x]=x;
 for(int a=1;a<A;++a)for(int x=0;x<X;++x)phi_table[a*X+x]=phi_table[(a-1)*X+x]-phi_table[(a-1)*X+x/ps[a-1]];
 for(int i=1;i<argc;++i){u64 n=std::stoull(argv[i]);if(n>(1ull<<36))return 2;std::cout<<n<<' '<<count(n)<<'\n';}
}
