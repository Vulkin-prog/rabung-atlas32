#include "../src/periodic_bands.hpp"
#include "../baseline/src/anchor_geometry.hpp"
#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>
using namespace vdw_periodic;
using Key=std::tuple<u32,u32,u32>;
void check(bool b,const char*m){if(!b)throw std::runtime_error(m);}
int main(){try{u64 configs=0,compared=0,periods=0,empty=0,poison=0;
 for(u32 B:{8u,17u,31u,32u,33u,63u,64u,65u,127u,255u,513u,1025u})
 for(u32 q:{3u,4u,5u,8u,11u,16u,25u})if(q<B)
 for(u32 t=3;t<=q;++t){u32 s=q-t+1,D=2*B/(q-1);
  for(u32 d=1;d<=D;++d){u32 g=gcd(d,6);for(u32 c=g==1?0:1;c<g;++c)if(g==1||gcd(c,g)==1)
   for(u32 phase=0;phase<s;++phase){++configs;
    std::vector<Key>want,got;
    for(auto b:vdw_anchor::bands(B,d,t,s,phase)){if(b.last<c)continue;u32 f=b.first>c?(b.first-c+g-1)/g:0,l=(b.last-c)/g;
     if(f<=l)for(u32 w=f/32;w<=l/32;++w)want.emplace_back(f,l,w);}
    auto pp=prepare(make_shapes(B,d,g,c,t,s,phase));
    for(u32 i=0;i<pp.total;++i){auto x=select(pp,i);check(x.valid&&x.first<=x.last&&x.word>=x.first/32&&x.word<=x.last/32,"invalid ranked task");got.emplace_back(x.first,x.last,x.word);++compared;}
    check(!select(pp,pp.total).valid,"out-of-range select accepted");
    std::sort(want.begin(),want.end());std::sort(got.begin(),got.end());check(want==got,"rank/select differs from enumerated bands");empty+=got.empty();
    for(auto pr:{pp.a,pp.b})if(pr.shape.middle){for(u32 k=0;k<pr.P;++k){check(band_words(pr.shape,k)==band_words(pr.shape,k+pr.P),"period not exact");++periods;}
     // An invalid assumed period=1 MUST be detected when lengths vary.
     bool varies=false;for(u32 k=1;k<pr.P;++k)varies|=band_words(pr.shape,k)!=band_words(pr.shape,0);if(varies)++poison;
    }
   }
  }
 }
 check(poison>0,"no nonconstant-period cases exercised");
 std::cout<<"{\"status\":\"PASS\",\"configurations\":"<<configs<<",\"tasks_compared\":"<<compared<<",\"period_equalities\":"<<periods<<",\"empty_plans\":"<<empty<<",\"constant_period_negative_cases\":"<<poison<<",\"GPU_executed\":false}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
