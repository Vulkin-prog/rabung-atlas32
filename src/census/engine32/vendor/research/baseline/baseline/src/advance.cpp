#define VDW_PERIODIC_NO_MAIN
#include "../previous/src/periodic_reference.cpp"
#include <map>
#include "coupled_device.hpp"

struct Extra {u64 words=0,stages=0,raw=0,lanes=0,rows=0;std::map<u32,u64> gate_words;std::vector<std::array<u32,3>> improvements;};
static void certify_witness(u32 p,u32 r,const AnchorStats& R,bool maximal){
 require(R.maximum>0&&R.witness>0&&(u64)R.witness+R.maximum<=p,"witness range");
 auto colour=power(R.witness,(p-1)/r,p);
 for(u32 j=0;j<R.maximum;++j)require(power(R.witness+j,(p-1)/r,p)==colour,"witness direct exponentiation");
 if(maximal&&R.witness>1)require(power(R.witness-1,(p-1)/r,p)!=colour,"left neighbour differs");
 if(maximal&&(u64)R.witness+R.maximum<p)require(power(R.witness+R.maximum,(p-1)/r,p)!=colour,"right neighbour differs");
}
struct Prepared {Table pilot_table;AnchorStats pilot;u32 t;
 Prepared(u32 p,u32 r,u32 pd):pilot_table(p,r,std::max(3u,select_t(p,r).t)),t(pilot_table.L){
  seed_table(pilot_table,pilot);if(pd)pilot=scan_bands(pilot_table,t,t,std::min(pd,pilot_table.D),false,false,true,pilot);
  require(pilot.maximum>=t,"fallback required");certify_witness(p,r,pilot,false);
 }
};
struct Output {AnchorStats R;Extra X;u32 t,Q,cap,B,D;u64 pilot_words;double seconds;};

template<class Visitor>
static Extra traverse(const Table&T,u32 t,u32 fixedQ,AnchorStats&R,bool adaptive,Visitor visit,u32 batch=1){
 Positive P(T);Wheel W2(T,2),W3(T,3),W6(T,6);auto inv=inverses(T.p,T.D);Extra X;
 u32 gate=fixedQ;require(batch>0,"positive epoch size");
 for(u32 d=1;d<=T.D;++d){
  // Snapshot once per denominator: conjugate residue rows share this gate.
  if((d-1)%batch==0)gate=adaptive?std::min(T.L,R.maximum):fixedQ;require(gate>=t&&gate<=T.L,"gate domain");
  u32 spacing=gate-t+1,g=std::gcd(d,6u);
  auto ph=vdw_anchor::choose(T.B,d,t,spacing,true,true);
  auto row=[&](auto view,u32 c){auto plan=vdw_periodic::prepare(vdw_periodic::make_shapes(T.B,d,g,c,t,spacing,ph.phase));
   ++X.rows;X.lanes+=32ull*((plan.total+31)/32);X.gate_words[gate]+=plan.total;
   for(u32 rank=0;rank<plan.total;++rank){auto task=vdw_periodic::select(plan,rank);++X.words;vdw_next::Counters st;
    u32 bits=vdw_next::filter_range(view,d/g,t,task.word,task.first,task.last,&st);X.stages+=st.stages;X.raw+=__builtin_popcount(bits);
    while(bits){u32 j=__builtin_ctz(bits);bits&=bits-1;u32 old=R.maximum;visit(g*(task.word*32+j)+c,d,inv[d]);if(R.maximum>old)X.improvements.push_back({d,R.maximum,R.witness});}
   }
  };
  if(g==1)row(P.view(),0);else{auto&W=g==2?W2:g==3?W3:W6;for(u32 k=0;k<W.residues.size();++k)row(W.view(k),W.residues[k]);}
 }
 return X;
}
static Output adaptive_run(u32 p,u32 r,u32 pd,u32 surplus,bool adaptive){require(surplus<=256,"surplus domain");auto begin=std::chrono::steady_clock::now();Prepared A(p,r,pd);u32 q=A.pilot.maximum,cap=q+surplus;
 require(cap<p&&cap<=256,"coverage cap domain");Table T(p,r,cap);extend_table(T,vdw_clipped::from_denominator(p,cap,std::max(1u,T.D/2)));
 AnchorStats R;R.maximum=A.pilot.maximum;R.witness=A.pilot.witness;seed_table(T,R);
 auto X=traverse(T,A.t,q,R,adaptive,[&](u32 b,u32 d,u32 inv){complete_anchor(T,A.t,b,d,inv,R);});
 if(p<=2000000){NewResult N;N.result.M=R.maximum;N.result.prefix=R.prefix;N.result.witness=R.witness;compare_dense_profile(T,N);}
 certify_witness(p,r,R,true);
 return {R,X,A.t,q,cap,T.B,T.D,A.pilot.words,std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()};
}
struct Run {u32 first,length; i64 numerator; u32 denominator;};
static bool full_run(const Table&T,u32 t,u32 b,u32 d,u32 inv,Run&run,u64&extensions){
 if(!vdw_next::admissible(b,d,t))return false;
 u32 x=(T.p-(u64)b*inv%T.p)%T.p;require(x>0&&(u64)x+t<=T.p,"lift");
 int col=(T.colour(-(i64)b)+T.r-T.c[d])%T.r;Stats st;
 auto at=[&](u32 xx)->int{i64 z=-(i64)b+((i64)xx-x)*d;if(z&&z>=-(i64)T.B&&z<=T.B)return(T.colour(z)+T.r-T.c[d])%T.r;return small_oracle(T,xx,st);};
 u32 l=x,h=x+t-1;while(l>1){++extensions;if(at(l-1)!=col)break;--l;}while(h<T.p-1){++extensions;if(at(h+1)!=col)break;++h;}
 i64 a=-(i64)b+((i64)l-x)*d;u32 len=h-l+1;
 if(T.p-h<l){l=T.p-h;a=-a-(i64)(len-1)*d;}
 run={l,len,a,d};return true;
}
static void print(const Output&O){std::cout<<"{\"status\":\"PASS\",\"t\":"<<O.t<<",\"Q\":"<<O.Q<<",\"cap\":"<<O.cap<<",\"B\":"<<O.B<<",\"D\":"<<O.D<<",\"maximum\":"<<O.R.maximum<<",\"prefix\":"<<O.R.prefix<<",\"witness\":"<<O.R.witness<<",\"words\":"<<O.X.words<<",\"pilot_words\":"<<O.pilot_words<<",\"stages\":"<<O.X.stages<<",\"raw\":"<<O.X.raw<<",\"slots\":"<<O.X.lanes<<",\"extensions\":"<<O.R.extensions<<",\"gate_words\":{";bool first=true;for(auto [g,w]:O.X.gate_words){if(!first)std::cout<<',';first=false;std::cout<<'"'<<g<<"\":"<<w;}std::cout<<"},\"improvements\":[";first=true;for(auto a:O.X.improvements){if(!first)std::cout<<',';first=false;std::cout<<'['<<a[0]<<','<<a[1]<<','<<a[2]<<']';}std::cout<<"],\"seconds\":"<<O.seconds<<"}\n";}
static void atlas_run(u32 p,u32 coarse,u32 fine,u32 t,u32 pd){auto begin=std::chrono::steady_clock::now();require(fine%coarse==0,"divisor");Prepared F(p,fine,pd);u32 K=F.pilot.maximum;require(K>=t,"threshold too small");
 Table T(p,coarse,K);extend_table(T,vdw_clipped::from_denominator(p,K,std::max(1u,T.D/2)));AnchorStats R;seed_table(T,R);std::map<u32,u32> atlas;u64 ext=0;
 auto X=traverse(T,t,K,R,false,[&](u32 b,u32 d,u32 inv){Run z;if(full_run(T,t,b,d,inv,z,ext)){update(R,p,z.first,z.length);if(z.length>=K){auto [it,inserted]=atlas.emplace(z.first,z.length);require(it->second==z.length,"atlas consistency");}}});
 Table Fine(p,fine,K);AnchorStats refined;refined.maximum=F.pilot.maximum;refined.witness=F.pilot.witness;seed_table(Fine,refined);Stats queries;u64 examined=0;
 for(auto [l,len]:atlas){u32 prev=UINT32_MAX,run=0;for(u32 j=0;j<len;++j){u32 c=small_oracle(Fine,l+j,queries);run=c==prev?run+1:1;prev=c;update(refined,p,l+j-run+1,run);++examined;}}
 if(p<=2000000){NewResult N;N.result.M=refined.maximum;N.result.prefix=refined.prefix;N.result.witness=refined.witness;compare_dense_profile(Fine,N);}
 std::cout<<"{\"status\":\"PASS\",\"p\":"<<p<<",\"coarse\":"<<coarse<<",\"fine\":"<<fine<<",\"t\":"<<t<<",\"K\":"<<K<<",\"B\":"<<T.B<<",\"D\":"<<T.D<<",\"coarse_maximum\":"<<R.maximum<<",\"fine_maximum\":"<<refined.maximum<<",\"fine_prefix\":"<<refined.prefix<<",\"fine_witness\":"<<refined.witness<<",\"words\":"<<X.words<<",\"stages\":"<<X.stages<<",\"raw\":"<<X.raw<<",\"extension_tests\":"<<ext<<",\"intervals\":"<<atlas.size()<<",\"refinement_positions\":"<<examined<<",\"oracle_euclid\":"<<queries.euclid<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<"}\n";
}

struct SharedOutput {u32 p,coarse,fine,t,K,cap,B,D;AnchorStats coarse_result,fine_result;Extra X;u64 pilot_words,refine_runs,refine_positions,outside,coarse_extensions;double seconds;};
static SharedOutput shared_run(u32 p,u32 coarse,u32 fine,u32 t,u32 pd,u32 surplus,u32 batch=1,bool check_core=false){auto begin=std::chrono::steady_clock::now();
 require(coarse>=2&&fine>coarse&&fine<=9&&fine%coarse==0&&surplus<=256,"proper refinement");Prepared F(p,fine,pd);u32 K=F.pilot.maximum,cap=K+surplus;if(!t)t=std::max(3u,K>=(coarse==2?4u:3u)?K-(coarse==2?4u:3u):3u);require(K>=t&&t>=3&&cap<=256&&cap<p,"shared domain");
 Table Fine(p,fine,cap);extend_table(Fine,vdw_clipped::from_denominator(p,cap,std::max(1u,Fine.D/2)));
 Table T=Fine;T.r=coarse;T.sigma%=coarse;for(auto& c:T.c)c%=coarse;
 u32 np=0;while((1u<<np)<coarse)++np;T.planes.assign(np,std::vector<u32>((T.n+31)/32+2));
 for(u32 i=0;i<T.n;++i)if(i!=T.B){int c=T.colour((i64)i-T.B);for(u32 j=0;j<np;++j)if((c>>j)&1)T.planes[j][i>>5]|=1u<<(i&31);}
 AnchorStats RC,RF;RF.maximum=RC.maximum=K;RF.witness=RC.witness=F.pilot.witness;seed_table(T,RC);seed_table(Fine,RF);
 u64 ext=0,runs=0,positions=0,outside=0;auto invers=inverses(p,T.D);vdw_next::CompletionView CV{p,fine,Fine.B,t,Fine.D,Fine.c.data(),invers.data()};
 auto X=traverse(T,t,K,RF,true,[&](u32 b,u32 d,u32 inv){
  vdw_pair::PairResult core{};
  if(check_core)core=vdw_pair::complete(CV,coarse,{d,b,1u,0u},vdw_next::pack_best(RC.maximum,RC.witness),vdw_next::pack_best(RF.maximum,RF.witness));
  auto verify=[&](){if(check_core)require(core.ok&&core.coarse_best==vdw_next::pack_best(RC.maximum,RC.witness)&&core.fine_best==vdw_next::pack_best(RF.maximum,RF.witness),"device-core/reference mismatch");};
  Run z;if(!full_run(T,t,b,d,inv,z,ext)){verify();return;}
  update(RC,p,z.first,z.length);if(z.length<RF.maximum){verify();return;}++runs;
  u32 prev=UINT32_MAX,run=0;Stats qs;
  for(u32 j=0;j<z.length;++j){i64 num=z.numerator+(i64)j*d;u32 c;
   if(num&&num>=-(i64)Fine.B&&num<=(i64)Fine.B)c=(Fine.colour(num)+fine-Fine.c[d])%fine;
   else{++outside;c=small_oracle(Fine,z.first+j,qs);}
   run=c==prev?run+1:1;prev=c;update(RF,p,z.first+j-run+1,run);++positions;
  }
  RF.queries+=qs.queries;RF.euclid+=qs.euclid;require(RC.maximum>=RF.maximum,"refinement invariant");verify();
 },batch);
 if(p<=2000000){for(auto pair:{std::make_pair(&T,&RC),std::make_pair(&Fine,&RF)}){NewResult N;N.result.M=pair.second->maximum;N.result.prefix=pair.second->prefix;N.result.witness=pair.second->witness;compare_dense_profile(*pair.first,N);}}
 certify_witness(p,coarse,RC,true);certify_witness(p,fine,RF,true);
 return {p,coarse,fine,t,K,cap,T.B,T.D,RC,RF,X,F.pilot.words,runs,positions,outside,ext,std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()};
}
static void print_shared(const SharedOutput&O){std::cout<<"{\"status\":\"PASS\",\"p\":"<<O.p<<",\"coarse\":"<<O.coarse<<",\"fine\":"<<O.fine<<",\"t\":"<<O.t<<",\"K\":"<<O.K<<",\"cap\":"<<O.cap<<",\"B\":"<<O.B<<",\"D\":"<<O.D
 <<",\"coarse_maximum\":"<<O.coarse_result.maximum<<",\"coarse_prefix\":"<<O.coarse_result.prefix<<",\"coarse_witness\":"<<O.coarse_result.witness<<",\"fine_maximum\":"<<O.fine_result.maximum<<",\"fine_prefix\":"<<O.fine_result.prefix<<",\"fine_witness\":"<<O.fine_result.witness
 <<",\"pilot_words\":"<<O.pilot_words<<",\"words\":"<<O.X.words<<",\"stages\":"<<O.X.stages<<",\"raw\":"<<O.X.raw<<",\"slots\":"<<O.X.lanes<<",\"coarse_extensions\":"<<O.coarse_extensions<<",\"refine_runs\":"<<O.refine_runs<<",\"refine_positions\":"<<O.refine_positions<<",\"outside_table\":"<<O.outside<<",\"gate_words\":{";bool first=true;for(auto[g,w]:O.X.gate_words){if(!first)std::cout<<',';first=false;std::cout<<'"'<<g<<"\":"<<w;}std::cout<<"},\"seconds\":"<<O.seconds<<"}\n";}
#ifndef VDW_ADVANCE_NO_MAIN
int main(int argc,char**argv){try{
 if(argc==6&&std::string(argv[1])=="separate"){
 auto p=parse(argv[2]),c=parse(argv[3]),f=parse(argv[4]),pd=parse(argv[5]);
 auto start=std::chrono::steady_clock::now();auto A=adaptive_run(p,c,pd,0,false),B=adaptive_run(p,f,pd,0,false);
 std::cout<<"{\"status\":\"PASS\",\"coarse_maximum\":"<<A.R.maximum<<",\"fine_maximum\":"<<B.R.maximum<<",\"coarse_prefix\":"<<A.R.prefix<<",\"fine_prefix\":"<<B.R.prefix<<",\"coarse_witness\":"<<A.R.witness<<",\"fine_witness\":"<<B.R.witness<<",\"words\":"<<A.X.words+B.X.words<<",\"pilot_words\":"<<A.pilot_words+B.pilot_words<<",\"stages\":"<<A.X.stages+B.X.stages<<",\"shifted_plane_reads\":"<<A.X.stages*(c<=2?1:c<=4?2:c<=8?3:4)+B.X.stages*(f<=2?1:f<=4?2:f<=8?3:4)<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";return 0;}
 if(argc==7&&std::string(argv[1])=="adaptive"){print(adaptive_run(parse(argv[2]),parse(argv[3]),parse(argv[4]),parse(argv[5]),parse(argv[6])));return 0;}
 if((argc==8||argc==9||argc==10)&&std::string(argv[1])=="shared"){print_shared(shared_run(parse(argv[2]),parse(argv[3]),parse(argv[4]),parse(argv[5]),parse(argv[6]),parse(argv[7]),argc>=9?parse(argv[8]):1,argc==10?parse(argv[9]):false));return 0;}
 if(argc==7&&std::string(argv[1])=="atlas"){atlas_run(parse(argv[2]),parse(argv[3]),parse(argv[4]),parse(argv[5]),parse(argv[6]));return 0;}
 throw std::invalid_argument("adaptive p r pilot surplus adaptive | atlas p coarse fine t pilot | shared p coarse fine t pilot surplus");
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}

#endif
