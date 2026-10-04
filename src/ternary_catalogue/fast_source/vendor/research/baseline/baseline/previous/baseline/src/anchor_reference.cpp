// Read-only reuse of the prior character table and word32 filter.
#define VDW_POSITIVE_NO_MAIN
#include "../vendor/src/positive_reference.cpp"
#include "anchor_geometry.hpp"
#include <iomanip>
struct AnchorStats {
 u64 words=0,stages=0,direct=0,signed_reads=0,raw=0,primitives=0,extensions=0,queries=0,euclid=0;
 u64 bands=0,phase_trials=0;u32 maximum=0,prefix=0,witness=0;
};
static void update(AnchorStats& R,u32 p,u32 start,u32 len){
 const u32 reflected=p-start-len+1;start=std::min(start,reflected);
 if(len>R.maximum||(len==R.maximum&&(!R.witness||start<R.witness))){R.maximum=len;R.witness=start;}
}
static void seed_table(const Table&T,AnchorStats&R){
 u32 run=0;int prev=-1;
 for(u32 x=1;x<=T.B;++x){int c=T.c[x];run=c==prev?run+1:1;prev=c;update(R,T.p,x-run+1,run);}
 Stats st;while(R.prefix+1<T.p&&small_oracle(T,R.prefix+1,st)==0)++R.prefix;
 update(R,T.p,1,R.prefix);R.queries+=st.queries;R.euclid+=st.euclid;
}
// IMPORTANT: these anchors need not be run HEADS. Both directions are completed.
// The reflected run has the same length, so it need not be completed a second time.
static void complete_anchor(const Table&T,u32 anchor_t,u32 b,u32 d,u32 inv,AnchorStats&R){
 if(!vdw_next::admissible(b,d,anchor_t))return;++R.primitives;
 const u32 x=(T.p-(u64)b*inv%T.p)%T.p;
 require(x>=1&&(u64)x+anchor_t<=T.p,"sampled anchor linear lift");
 const int color=(T.colour(-(i64)b)+T.r-T.c[d])%T.r;
 Stats st;
 auto at=[&](u32 xx)->int{
  const i64 z=-(i64)b+((i64)xx-x)*d;
  if(z&&z>=-(i64)T.B&&z<=(i64)T.B)return (T.colour(z)+T.r-T.c[d])%T.r;
  return small_oracle(T,xx,st);
 };
 u32 left=x,right=x+anchor_t-1;
 while(left>1){++R.extensions;if(at(left-1)!=color)break;--left;}
 while(right<T.p-1){++R.extensions;if(at(right+1)!=color)break;++right;}
 update(R,T.p,left,right-left+1);R.queries+=st.queries;R.euclid+=st.euclid;
}
static std::vector<u32> inverses(u32 p,u32 D){std::vector<u32> inv(D+1);inv[1]=1;for(u32 d=2;d<=D;++d)inv[d]=p-(u64)(p/d)*inv[p%d]%p;return inv;}
static AnchorStats scan_bands(const Table&T,u32 anchor_t,u32 coverage_q,u32 last_d,bool sample,bool optimize,bool wheel,const AnchorStats&seed){
 AnchorStats R=seed;Positive P(T);Wheel W2(T,2),W3(T,3),W6(T,6);auto inv=inverses(T.p,last_d);
 const u32 spacing=sample?coverage_q-anchor_t+1:1;
 for(u32 d=1;d<=last_d;++d){
  std::vector<vdw_anchor::Band> bands;
  if(sample){auto plan=vdw_anchor::choose(T.B,d,anchor_t,spacing,wheel,optimize);R.phase_trials+=optimize?spacing:1;bands=vdw_anchor::bands(T.B,d,anchor_t,spacing,plan.phase);}
  else{u32 first=((anchor_t-1)*d+1)/2;if(first<=T.B)bands.push_back({first,T.B});}
  R.bands+=bands.size();u32 g=wheel?std::gcd(d,6u):1;
  auto traverse=[&](auto view,u32 c,const vdw_anchor::Band&band){
   if(band.last<c)return;u32 first=band.first>c?(band.first-c+g-1)/g:0,maximum=(band.last-c)/g;
   if(first>maximum)return;
   for(u32 w=first>>5;w<=maximum/32;++w){
    ++R.words;vdw_next::Counters fs;
    u32 bits=vdw_next::filter_range(view,d/g,anchor_t,w,first,maximum,&fs);
    R.stages+=fs.stages;R.direct+=fs.direct;R.signed_reads+=fs.signed_reads;
    R.raw+=__builtin_popcount(bits);
    while(bits){u32 k=__builtin_ctz(bits);bits&=bits-1;complete_anchor(T,anchor_t,g*(32*w+k)+c,d,inv[d],R);}
   }
  };
  for(auto band:bands){
   if(g==1)traverse(P.view(),0,band);
   else{const Wheel&W=g==2?W2:g==3?W3:W6;for(u32 row=0;row<W.residues.size();++row)traverse(W.view(row),W.residues[row],band);}
  }
 }
 return R;
}
struct Experiment{u32 p,r,t,pilot_d,q,B,D;AnchorStats pilot,main;bool dense_checked=false;double seconds=0,main_seconds=0;};
static Experiment experiment(u32 p,u32 r,u32 pilot_d,bool optimize=true){
 require(prime(p)&&p>=7&&r>=2&&r<=9&&(p-1)%r==0,"experiment domain");
 auto begin=std::chrono::steady_clock::now();const u32 t=std::max(3u,select_t(p,r).t);
 Table P(p,r,t);AnchorStats initial;seed_table(P,initial);
 AnchorStats pil=initial;
 if(pilot_d)pil=scan_bands(P,t,t,std::min(pilot_d,P.D),false,false,true,initial);
 u32 q=pil.maximum;
 const auto pilot_colour=power(pil.witness,(p-1)/r,p);
 for(u32 j=0;j<q;++j)require(power(pil.witness+j,(p-1)/r,p)==pilot_colour,"pilot witness by direct exponentiation");
 require(q>=t,"pilot too short: dense fallback is required");
 Table T(p,r,q); // This table covers FULL q-term windows; the filter tests t terms.
 AnchorStats starting;starting.maximum=pil.maximum;starting.witness=pil.witness;seed_table(T,starting);
 const auto main_start=std::chrono::steady_clock::now();
 auto R=scan_bands(T,t,q,T.D,true,optimize,true,starting);
 const double main_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-main_start).count();
 bool checked=p<=2000000;
 if(checked){NewResult N;N.result.M=R.maximum;N.result.prefix=R.prefix;N.result.witness=R.witness;compare_dense_profile(T,N);}
 auto ch=power(R.witness,(p-1)/r,p);
 for(u32 j=0;j<R.maximum;++j)require(power(R.witness+j,(p-1)/r,p)==ch,"direct witness");
 if(R.witness>1)require(power(R.witness-1,(p-1)/r,p)!=ch,"left neighbour");
 if((u64)R.witness+R.maximum<p)require(power(R.witness+R.maximum,(p-1)/r,p)!=ch,"right neighbour");
 return {p,r,t,std::min(pilot_d,P.D),q,T.B,T.D,pil,R,checked,std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count(),main_seconds};
}
static void print_experiment(const Experiment&E){auto&R=E.main;auto&P=E.pilot;
 std::cout<<std::setprecision(8)<<"{\"status\":\"PASS\",\"p\":"<<E.p<<",\"r\":"<<E.r<<",\"anchor_t\":"<<E.t<<",\"pilot_last_d\":"<<E.pilot_d<<",\"certified_seed\":"<<E.q<<",\"seed_witness\":"<<P.witness<<",\"seed_checked_directly\":true,\"spacing\":"<<E.q-E.t+1<<",\"B\":"<<E.B<<",\"D\":"<<E.D
 <<",\"maximum\":"<<R.maximum<<",\"prefix\":"<<R.prefix<<",\"witness\":"<<R.witness
 <<",\"pilot_words\":"<<P.words<<",\"main_words\":"<<R.words<<",\"total_words\":"<<P.words+R.words<<",\"pilot_stages\":"<<P.stages<<",\"main_stages\":"<<R.stages<<",\"total_stages\":"<<P.stages+R.stages
 <<",\"main_bands\":"<<R.bands<<",\"phase_trials\":"<<R.phase_trials<<",\"raw_bits\":"<<R.raw<<",\"primitive_anchors\":"<<R.primitives<<",\"extensions\":"<<R.extensions<<",\"oracle_queries\":"<<R.queries<<",\"euclid_steps\":"<<R.euclid
 <<",\"main_plane_queries\":"<<R.direct+R.signed_reads<<",\"direct_plane_queries\":"<<R.direct
 <<",\"dense_checked\":"<<(E.dense_checked?"true":"false")<<",\"cpu_seconds_including_tables_and_planner\":"<<E.seconds<<",\"main_cpu_seconds_including_planner\":"<<E.main_seconds<<",\"GPU_executed\":false}\n";
}
#ifndef VDW_ANCHOR_NO_MAIN
int main(int argc,char**argv){try{if(argc!=5){std::cerr<<"usage: --profile PRIME R PILOT_D\n";return 2;}require(std::string(argv[1])=="--profile","mode");print_experiment(experiment(parse(argv[2]),parse(argv[3]),parse(argv[4])));return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
