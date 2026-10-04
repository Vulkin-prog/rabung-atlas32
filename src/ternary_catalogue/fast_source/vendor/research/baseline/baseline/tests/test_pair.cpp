#define VDW_ADVANCE_NO_MAIN
#include "../src/advance.cpp"
struct WordProfile{u32 maximum=0,start=0;};
static WordProfile plain(const std::vector<unsigned char>&w,u32 r){WordProfile z;u32 run=0,prev=100;for(u32 i=0;i<w.size();++i){u32 c=w[i]%r;run=c==prev?run+1:1;prev=c;if(run>z.maximum){z.maximum=run;z.start=i+1-run;}}return z;}
static u64 word_test(){u64 cases=0;
 for(auto [coarse,fine,limit]:{std::tuple{2u,4u,8u},std::tuple{3u,6u,6u},std::tuple{2u,8u,5u},std::tuple{3u,9u,5u}}){u64 total=1;
  for(u32 n=1;n<=limit;++n){total*=fine;std::vector<unsigned char>w(n);
   for(u64 rank=0;rank<total;++rank){u64 code=rank;for(auto&v:w){v=code%fine;code/=fine;}auto want=plain(w,fine);
    for(u32 K=1;K<=want.maximum;++K){u32 got=0;
     for(u32 l=0;l<n;){u32 h=l+1;while(h<n&&w[h]%coarse==w[l]%coarse)++h;
      if(h-l>=K){std::vector<unsigned char>part(w.begin()+l,w.begin()+h);got=std::max(got,plain(part,fine).maximum);}l=h;}
     require(got==want.maximum,"word refinement lemma");++cases;
    }
   }
  }
 }
 return cases;
}
static void negatives(){
 {u32 p=113;Table F(p,4,4);std::vector<unsigned char>w(p-1);u32 x=1;for(u32 j=0;j<p-1;++j){w[x-1]=j%4;x=(u64)x*3%p;}
 require(plain(w,2).maximum==5&&plain(w,4).maximum==4,"p113 full profiles");u32 restricted=0;
 for(u32 l=0;l<w.size();){u32 h=l+1;while(h<w.size()&&w[h]%2==w[l]%2)++h;
 if(h-l==5){std::vector<unsigned char>a(w.begin()+l,w.begin()+h);restricted=std::max(restricted,plain(a,4).maximum);}l=h;}
 require(restricted==3,"wrong coarse-gate counterexample");}
 Table T(1009,4,5);auto inv=inverses(T.p,T.D);vdw_next::CompletionView V{T.p,T.r,T.B,3,T.D,T.c.data(),inv.data()};u64 seed=vdw_next::pack_best(1,1);
 require(!vdw_pair::complete(V,3,{1,4,1,0},seed,seed).ok,"non-divisor rejection");
 require(!vdw_pair::complete(V,2,{0,4,1,0},seed,seed).ok,"zero denominator rejection");
 require(!vdw_pair::complete(V,2,{1,4,1,1},seed,seed).ok,"incompatible descriptor rejection");
 require(!vdw_pair::complete(V,2,{1,T.B+1,1,0},seed,seed).ok,"out of range rejection");
 V.colour=nullptr;require(!vdw_pair::complete(V,2,{1,4,1,0},seed,seed).ok,"null table rejection");
}
int main(int argc,char**argv){try{std::string mode=argc>1?argv[1]:"small";u64 counts=0,fallback=0;
 if(mode=="words"){std::cout<<"{\"status\":\"PASS\",\"word_threshold_cases\":"<<word_test()<<"}\n";return 0;}
 if(mode=="negative"){negatives();std::cout<<"{\"status\":\"PASS\",\"negative_cases\":6}\n";return 0;}
 std::vector<u32> ps;if(mode=="small"){for(u32 p=101;p<=1600;++p)if(prime(p))ps.push_back(p);}else if(mode=="million"){for(u32 p=1000000;p<1000600&&ps.size()<16;++p)if(prime(p))ps.push_back(p);}else throw std::invalid_argument("test mode");
 for(auto p:ps)for(auto pair:{std::pair{2u,4u},std::pair{2u,8u},std::pair{4u,8u},std::pair{3u,6u},std::pair{3u,9u}})if((p-1)%pair.second==0){
  // Explicitly recognize only the certified lack of a sufficiently long seed.
  u32 q=0;try{Prepared F(p,pair.second,8);q=F.pilot.maximum;}catch(const std::exception&e){if(std::string(e.what())=="fallback required"){++fallback;continue;}throw;}
  for(u32 surplus:{0u,3u})for(u32 batch:{1u,64u}){
   auto O=shared_run(p,pair.first,pair.second,std::min(3u,q),8,surplus,batch,true);require(O.coarse_result.maximum>=O.fine_result.maximum,"nested maxima");++counts;
  }
 }
 std::cout<<"{\"status\":\"PASS\",\"mode\":\""<<mode<<"\",\"coupled_runs_checked_against_dense\":"<<counts<<",\"individual_profile_checks\":"<<2*counts<<",\"seed_fallbacks\":"<<fallback<<",\"device_core_compared_per_anchor\":true}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
