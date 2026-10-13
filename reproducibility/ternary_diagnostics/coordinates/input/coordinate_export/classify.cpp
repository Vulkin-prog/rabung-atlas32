// Read-only classification of retained ordinary runs; no character evaluation.
// Domain: p > 8192, 1 <= x <= (p-1)/2, x+length-1 <= (p-1)/2.
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
using U=uint64_t;
struct QHit {int j=-1,A=0,D=0;};
static bool nonsquare(int d) {
 if(d<=0) return d<0;
 for(int t=1;t*t<=d;++t) if(t*t==d) return false;
 return true;
}
struct Hit {QHit small,wide; int q=0,a=0;};
static Hit classify(U p,U x,int len) {
 if(p<=8192 || p >= (U(1)<<34) || len<1 || len>65535 || x<1 || x>(p-1)/2 || U(len-1)>(p-1)/2-x) throw std::runtime_error("invalid run domain");
 Hit h; static const auto allowed=[](){std::array<bool,8193> a{};for(int d=-4096;d<=4096;++d)a[d+4096]=nonsquare(d);return a;}();
 int64_t z0=(int64_t)(2*x)-128; if(z0<0) z0+=(int64_t)p;
 U z=(U)z0, v=(U)((__uint128_t)z*z%p), step=(2*z+1)%p;
 for(int b=-128;b<=128+2*(len-1);++b) {
  bool near=(v<=4096 || v>=p-4096);
  if(near){
   int d=(v<=4096)?(int)v:(int)((int64_t)v-(int64_t)p);
   if(allowed[d+4096] && (((int64_t)b*b-d)%4==0)) {
    if(h.wide.j<0){int j=b>128?(b-128+1)/2:0;h.wide={j,b-2*j,d};}
    if(h.small.j<0 && d>=-256 && d<=256 && b>=-64 && b<=64+2*(len-1)){
      int j=b>64?(b-64+1)/2:0;h.small={j,b-2*j,d};
    }
   }
  }
  v+=step; if(v>=p) v-=p; step+=2; if(step>=p) step-=p;
 }
 U residue=x;
 for(int q=2;q<=128;++q){residue+=x;if(residue>=p)residue-=p;
  if(residue+U(q)*U(len-1)>=p){int a=(int)((__uint128_t)q*x/p)+1; if(std::gcd(a,q)==1){h.q=q;h.a=a;break;}}
 }
 return h;
}
int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("usage: classify INPUT.tsv HITS.csv SUMMARY.tsv");
 std::ifstream in(argv[1]);std::ofstream out(argv[2]),sum(argv[3]);
 if(!in||!out||!sum)throw std::runtime_error("open failed");
 out<<"stratum,p,start,length,colour,mask,small_j,small_A,small_D,wide_j,wide_A,wide_D,q,a\n";
 // categories: small quadratic, wide quadratic, rational, wide union rational
 U total[4][7]={},removed[4][7][4]={},affected[4][7][4]={},covered[4][7][4]={};
 U currentp=0, orig[7]={}, rem[7][4]={};int currents=-1;U rows=0;
 auto finish=[&](){if(currents<0)return;for(int t=0;t<7;++t)for(int c=0;c<4;++c){if(rem[t][c])++affected[currents][t][c];if(orig[t]&&orig[t]==rem[t][c])++covered[currents][t][c];}for(int t=0;t<7;++t){orig[t]=0;for(int c=0;c<4;++c)rem[t][c]=0;}};
 U p,x;int s,len,col;
 while(in>>s>>p>>x>>len>>col){
  if(s<0||s>3||col<0||col>2||len<21)throw std::runtime_error("invalid record");
  if(p!=currentp||s!=currents){if(currents>s||(currents==s&&p<currentp))throw std::runtime_error("unsorted input");finish();currentp=p;currents=s;}
  Hit h=classify(p,x,len); int mask=(h.small.j>=0?1:0)|(h.wide.j>=0?2:0)|(h.q?4:0);
  if((mask&1)&&!(mask&2))throw std::runtime_error("non-nested catalogues");
  bool tags[4]={bool(mask&1),bool(mask&2),bool(mask&4),bool(mask&6)};
  for(int t=0;t<7;++t)if(len>=21+t){++orig[t];++total[s][t];for(int c=0;c<4;++c)if(tags[c]){++removed[s][t][c];++rem[t][c];}}
  if(mask)out<<s<<','<<p<<','<<x<<','<<len<<','<<col<<','<<mask<<','<<h.small.j<<','<<h.small.A<<','<<h.small.D<<','<<h.wide.j<<','<<h.wide.A<<','<<h.wide.D<<','<<h.q<<','<<h.a<<'\n';
  ++rows;
 }
 if(!in.eof())throw std::runtime_error("malformed input");
 finish();
 sum<<"stratum,threshold,category,ordinary_runs,tagged_runs,affected_profiles,fully_covered_profiles\n";
 const char* names[4]={"quadratic_small","quadratic_wide","rational","union_wide"};
 for(int a=0;a<4;++a)for(int t=0;t<7;++t)for(int c=0;c<4;++c)sum<<a<<','<<21+t<<','<<names[c]<<','<<total[a][t]<<','<<removed[a][t][c]<<','<<affected[a][t][c]<<','<<covered[a][t][c]<<'\n';
 if(!out||!sum)throw std::runtime_error("write failed");
 std::cerr<<"classified ordinary runs: "<<rows<<"\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
