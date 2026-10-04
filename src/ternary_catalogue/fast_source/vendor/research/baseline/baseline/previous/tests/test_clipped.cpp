#include "../src/clipped_cover.hpp"
#include "../baseline/src/anchor_geometry.hpp"
#include <iostream>
#include <numeric>
#include <cmath>
using namespace vdw_clipped;
bool prime(u32 p){for(u32 q=2;q*q<=p;++q)if(p%q==0)return false;return p>=2;}
void check(bool b,const char*m){if(!b)throw std::runtime_error(m);}
int main(){try{u64 covers=0,windows=0,anchors=0,monotone=0,negative=0;
 for(u32 p=7;p<=251;p+=2)if(prime(p))for(u32 Q=3;Q<=19&&Q<p;++Q)
 for(u32 div:{1u,2u,3u}){
  u32 D=std::max(1u,(u32)(2*std::sqrt((double)p/(Q-1)))/div);Cover c;
  try{c=from_denominator(p,Q,D);}catch(const std::invalid_argument&){continue;}
  ++covers;auto min=from_table(p,Q,c.B);check(certified(p,min)&&min.D<=D,"minimum denominator");if(min.D>1)check(!certified(p,{c.B,min.D-1,Q}),"minimum D not minimal");++monotone;
  if(c.B>1&&!certified(p,{c.B-1,c.D,c.Q}))++negative;
  for(u32 x=1;x+Q<=p;++x){bool covered=false;int witness_a=0;u32 witness_d=0;
   for(u32 d=1;d<=c.D&&!covered;++d)for(int a=-(int)c.B;a+(int)((Q-1)*d)<=(int)c.B;++a){
    if(((std::int64_t)a-(std::int64_t)d*x)%p||std::gcd((u32)std::abs(a),d)!=1)continue;
    witness_a=a;witness_d=d;covered=true;break;
   }
   check(covered,"uncovered global window");++windows;
   for(u32 t=3;t<=Q;t+=2){u32 s=Q-t+1;for(u32 phase=0;phase<s;++phase){bool hit=false;
    for(u32 j=0;j<s;++j){std::int64_t b=-(std::int64_t)witness_a-(std::int64_t)j*witness_d;
     if(!vdw_anchor::selected(b,witness_d,t,s,phase))continue;
     if(2*b<(std::int64_t)(t-1)*witness_d)b=(std::int64_t)(t-1)*witness_d-b;
     check(b>=0&&b<=c.B&&vdw_anchor::selected(b,witness_d,t,s,phase),"anchor lost after clipped reflection");hit=true;
    }
    check(hit,"selected anchor missing");++anchors;
   }}
  }
 }
 check(!certified(101,{1,1,3}),"negative cover accepted");
 std::cout<<"{\"status\":\"PASS\",\"covers\":"<<covers<<",\"global_windows\":"<<windows<<",\"anchor_phase_checks\":"<<anchors<<",\"minimum_D_checks\":"<<monotone<<",\"B_minus_one_negative\":"<<negative<<",\"GPU_executed\":false}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
