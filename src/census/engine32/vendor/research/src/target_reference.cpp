#define VDW_CUT_NO_MAIN
#include "../baseline/src/cut_reference.cpp"
#include "profile_masks.hpp"
static u32 encoded(u32 colour,u32 r){u32 two=1;while(r%(2*two)==0)two*=2;return colour%two+two*(colour%(r/two));}
// Minimal number of binary coordinates separating each coarse fibre.
static u32 minimal_planes(u32 c,u32 f){
 require(c>=2&&f>c&&f<=9&&f%c==0,"refinement");u32 np=0;while((1u<<np)<f)++np;
 u32 best=(1u<<np)-1;
 for(u32 s=1;s<(1u<<np);++s){bool ok=true;
  for(u32 a=0;a<f;++a)for(u32 b=a+1;b<f;++b)if(a%c==b%c&&((encoded(a,f)^encoded(b,f))&s)==0)ok=false;
  if(ok&&(__builtin_popcount(s)<__builtin_popcount(best)||(__builtin_popcount(s)==__builtin_popcount(best)&&s<best)))best=s;
 }
 return best;
}
struct TargetCounts {vdw_target::Counters masks;u64 coarse_visits=0,fine_visits=0,fine_extensions=0;};
template<class Visitor>
static Extra traverse_targets(const Table&T,const Table&F,u32 t,u32 fixedQ,AnchorStats&RC,AnchorStats&RF,Visitor visit,u32 batch,TargetCounts&cnt){
 Positive P(T),FP(F);Wheel W2(T,2),W3(T,3),W6(T,6),FW2(F,2),FW3(F,3),FW6(F,6);auto inv=inverses(T.p,T.D);Extra X;
 const u32 selected=minimal_planes(T.r,F.r);u32 gate=fixedQ,qc=RC.maximum,qf=RF.maximum;require(batch>0,"epoch size");
 for(u32 d=1;d<=T.D;++d){
  if((d-1)%batch==0){qc=RC.maximum;qf=RF.maximum;gate=std::min({T.L,qc,qf});}require(gate>=t&&gate<=T.L,"gate");
  u32 spacing=gate-t+1,g=std::gcd(d,6u);auto ph=vdw_anchor::choose(T.B,d,t,spacing,true,true);
  auto row=[&](auto view,auto fineview,u32 c){auto plan=vdw_periodic::prepare(vdw_periodic::make_shapes(T.B,d,g,c,t,spacing,ph.phase));
   ++X.rows;X.lanes+=32ull*((plan.total+31)/32);X.gate_words[gate]+=plan.total;
   for(u32 rank=0;rank<plan.total;++rank){auto task=vdw_periodic::select(plan,rank);++X.words;vdw_next::Counters st;
    u32 bits=vdw_next::filter_range(view,d/g,t,task.word,task.first,task.last,&st);X.stages+=st.stages;X.raw+=__builtin_popcount(bits);
    auto mask=vdw_target::split(view,fineview,selected,d/g,t,qc,qf,task.word,bits,&cnt.masks);
    u32 rest=mask.coarse|mask.fine;
    while(rest){u32 bit=__builtin_ctz(rest);rest&=rest-1;bool coarse=(mask.coarse>>bit)&1;
     u32 old=RF.maximum;visit(g*(task.word*32+bit)+c,d,inv[d],coarse);
     if(RF.maximum>old)X.improvements.push_back({d,RF.maximum,RF.witness});
    }
   }
  };
  if(g==1)row(P.view(),FP.view(),0);else{auto&W=g==2?W2:g==3?W3:W6;auto&FW=g==2?FW2:g==3?FW3:FW6;for(u32 k=0;k<W.residues.size();++k)row(W.view(k),FW.view(k),W.residues[k]);}
 }
 return X;
}
static SharedOutput targeted_run(u32 p,u32 coarse,u32 fine,u32 t,u32 pd,u32 surplus,u32 batch,TargetCounts&cnt){
 auto begin=std::chrono::steady_clock::now();require(coarse>=2&&fine>coarse&&fine<=9&&fine%coarse==0&&surplus<=256,"proper refinement");Prepared F(p,fine,pd);u32 K=F.pilot.maximum,cap=K+surplus;
 if(t==1)t=std::max(3u,(K+1)/2);if(!t)t=std::max(3u,K>=(coarse==2?4u:3u)?K-(coarse==2?4u:3u):3u);require(K>=t&&t>=3&&cap<=256&&cap<p,"target domain");
 Table Fine(p,fine,cap);extend_table(Fine,vdw_clipped::from_denominator(p,cap,std::max(1u,Fine.D/2)));
 Table T=Fine;T.r=coarse;T.sigma%=coarse;for(auto&c:T.c)c%=coarse;
 u32 np=0;while((1u<<np)<coarse)++np;T.planes.assign(np,std::vector<u32>((T.n+31)/32+2));
 for(u32 i=0;i<T.n;++i)if(i!=T.B){int c=T.colour((i64)i-T.B);for(u32 j=0;j<np;++j)if((c>>j)&1)T.planes[j][i>>5]|=1u<<(i&31);}
 AnchorStats RC,RF;RF.maximum=RC.maximum=K;RF.witness=RC.witness=F.pilot.witness;seed_table(T,RC);seed_table(Fine,RF);
 u64 ext=0,runs=0,positions=0,outside=0;
 auto X=traverse_targets(T,Fine,t,K,RC,RF,[&](u32 b,u32 d,u32 inv,bool coarse_path){
  Run z;
  if(!coarse_path){
   if(!full_run(Fine,t,b,d,inv,z,cnt.fine_extensions))return;
   ++cnt.fine_visits;update(RF,p,z.first,z.length);require(RC.maximum>=RF.maximum,"fine-only upper invariant");return;
  }
  if(!full_run(T,t,b,d,inv,z,ext))return;
  ++cnt.coarse_visits;update(RC,p,z.first,z.length);if(z.length<RF.maximum)return;++runs;
  u32 prev=UINT32_MAX,run=0;Stats qs;
  for(u32 j=0;j<z.length;++j){i64 num=z.numerator+(i64)j*d;u32 c;
   if(num&&num>=-(i64)Fine.B&&num<=(i64)Fine.B)c=(Fine.colour(num)+fine-Fine.c[d])%fine;
   else{++outside;c=small_oracle(Fine,z.first+j,qs);}
   run=c==prev?run+1:1;prev=c;update(RF,p,z.first+j-run+1,run);++positions;
  }
  RF.queries+=qs.queries;RF.euclid+=qs.euclid;require(RC.maximum>=RF.maximum,"refinement invariant");
 },batch,cnt);
 if(p<=2000000){for(auto pair:{std::make_pair(&T,&RC),std::make_pair(&Fine,&RF)}){NewResult N;N.result.M=pair.second->maximum;N.result.prefix=pair.second->prefix;N.result.witness=pair.second->witness;compare_dense_profile(*pair.first,N);}}
 certify_witness(p,coarse,RC,true);certify_witness(p,fine,RF,true);
 return {p,coarse,fine,t,K,cap,T.B,T.D,RC,RF,X,F.pilot.words,runs,positions,outside,ext,std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()};
}
#ifndef VDW_TARGET_NO_MAIN
int main(int argc,char**argv){try{
 require(argc==8,"usage: target p coarse fine pilot mode batch t (mode 0=baseline, 1=targeted)");
 u32 p=parse(argv[1]),c=parse(argv[2]),f=parse(argv[3]),pd=parse(argv[4]),mode=parse(argv[5]),batch=parse(argv[6]),t=parse(argv[7]);require(mode<=1,"mode");
 TargetCounts cnt;vdw_cut::Counts cut;
 auto O=mode?targeted_run(p,c,f,t,pd,3,batch,cnt):shared_screened(p,c,f,t,pd,3,batch,false,1,&cut);
 std::ostringstream tmp;auto* old=std::cout.rdbuf(tmp.rdbuf());print_shared(O);std::cout.rdbuf(old);
 auto s=tmp.str();s.erase(s.rfind('}'));
 std::cout<<s<<",\"mode\":"<<mode<<",\"selected_planes\":"<<minimal_planes(c,f)<<",\"target_words\":"<<cnt.masks.words<<",\"coarse_bits\":"<<cnt.masks.coarse_bits<<",\"fine_bits\":"<<cnt.masks.fine_bits<<",\"residual_queries\":"<<cnt.masks.residual_queries<<",\"coarse_cut_queries\":"<<cnt.masks.coarse_screen_queries<<",\"fine_cut_queries\":"<<cnt.masks.fine_screen_queries<<",\"coarse_visits\":"<<cnt.coarse_visits<<",\"fine_visits\":"<<cnt.fine_visits<<",\"fine_extensions\":"<<cnt.fine_extensions<<",\"baseline_cut_output\":"<<cut.output<<",\"gpu_executed\":false}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
