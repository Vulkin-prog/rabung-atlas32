#include "../src/warp_select.hpp"
#include <array>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace vdw_periodic;
void check(bool b,const char*m){if(!b)throw std::runtime_error(m);}
int main(){try{u64 plans=0,checks=0,collective_vectors=0;std::mt19937 gen(47291);
 for(u32 i=0;i<20000;++i){u32 B=32+gen()%1048545,Q=3+gen()%254,t=3+gen()%(Q-2),s=Q-t+1,D=2*B/(Q-1);if(!D)continue;u32 d=1+gen()%D,g=gcd(d,6),c=g==1?0:((g>2&&gen()%2)?g-1:1),phase=gen()%s;
  auto sides=make_shapes(B,d,g,c,t,s,phase);
  for(auto shape:{sides.a,sides.b}){auto ref=prepare(shape);std::array<u32,32> prefix{};
   for(u32 lane=0;lane<32;++lane)prefix[lane]=lane<ref.P?(shape.middle?band_words(shape,lane):1):0;
   for(u32 delta=1;delta<32;delta*=2){auto old=prefix;for(u32 lane=delta;lane<32;++lane)prefix[lane]+=old[lane-delta];++collective_vectors;}
   for(u32 j=0;j<ref.P;++j)check(prefix[j]==ref.prefix[j+1],"warp prefix scan mismatch");
   DevicePlan p{shape,ref.P,prefix[ref.P-1],ref.mid_words,ref.total};++plans;
   for(u32 block=0;block<3;++block){u32 base=block==0?0:block==1?(ref.total?gen()%ref.total:0):(ref.total/32)*32;
    for(u32 lane=0;lane<32;++lane){u32 rank=base+lane;auto got=warp_select(p,rank,[&](u32 src){check(src<32,"shuffle invalid source");return prefix[src];});auto want=select(ref,rank);
     check(got.valid==want.valid&&got.first==want.first&&got.last==want.last&&got.word==want.word,"device rank decoder mismatch");++checks;}
   }
  }
 }
 std::cout<<"{\"status\":\"PASS\",\"plans\":"<<plans<<",\"decoded_lane_tasks\":"<<checks<<",\"prefix_scan_vectors\":"<<collective_vectors<<",\"GPU_executed\":false}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
