#define VDW_ADVANCE_NO_MAIN
#include "../baseline/src/advance.cpp"
#include "two_cuts.hpp"
template<class Visitor>
static Extra traverse_screened(const Table&T,u32 t,u32 fixedQ,AnchorStats&R,bool adaptive,Visitor visit,u32 batch, unsigned cut_mode, vdw_cut::Counts& cut){
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
    bits=vdw_cut::screen(view,d/g,t,gate,task.word,bits,cut_mode,&cut);
    while(bits){u32 j=__builtin_ctz(bits);bits&=bits-1;u32 old=R.maximum;visit(g*(task.word*32+j)+c,d,inv[d]);if(R.maximum>old)X.improvements.push_back({d,R.maximum,R.witness});}
   }
  };
  if(g==1)row(P.view(),0);else{auto&W=g==2?W2:g==3?W3:W6;for(u32 k=0;k<W.residues.size();++k)row(W.view(k),W.residues[k]);}
 }
 return X;
}
static SharedOutput shared_screened(u32 p,u32 coarse,u32 fine,u32 t,u32 pd,u32 surplus,u32 batch=1,bool check_core=false,unsigned cut_mode=1,vdw_cut::Counts* cut_stats=nullptr){
 vdw_cut::Counts local_cut; auto& cut=cut_stats?*cut_stats:local_cut;auto begin=std::chrono::steady_clock::now();
 require(coarse>=2&&fine>coarse&&fine<=9&&fine%coarse==0&&surplus<=256,"proper refinement");Prepared F(p,fine,pd);u32 K=F.pilot.maximum,cap=K+surplus;if(t==1)t=std::max(3u,(K+1)/2);if(!t)t=std::max(3u,K>=(coarse==2?4u:3u)?K-(coarse==2?4u:3u):3u);require(K>=t&&t>=3&&cap<=256&&cap<p,"shared domain");
 Table Fine(p,fine,cap);extend_table(Fine,vdw_clipped::from_denominator(p,cap,std::max(1u,Fine.D/2)));
 Table T=Fine;T.r=coarse;T.sigma%=coarse;for(auto& c:T.c)c%=coarse;
 u32 np=0;while((1u<<np)<coarse)++np;T.planes.assign(np,std::vector<u32>((T.n+31)/32+2));
 for(u32 i=0;i<T.n;++i)if(i!=T.B){int c=T.colour((i64)i-T.B);for(u32 j=0;j<np;++j)if((c>>j)&1)T.planes[j][i>>5]|=1u<<(i&31);}
 AnchorStats RC,RF;RF.maximum=RC.maximum=K;RF.witness=RC.witness=F.pilot.witness;seed_table(T,RC);seed_table(Fine,RF);
 u64 ext=0,runs=0,positions=0,outside=0;auto invers=inverses(p,T.D);vdw_next::CompletionView CV{p,fine,Fine.B,t,Fine.D,Fine.c.data(),invers.data()};
 auto X=traverse_screened(T,t,K,RF,true,[&](u32 b,u32 d,u32 inv){
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
 },batch,cut_mode,cut);
 if(p<=2000000){for(auto pair:{std::make_pair(&T,&RC),std::make_pair(&Fine,&RF)}){NewResult N;N.result.M=pair.second->maximum;N.result.prefix=pair.second->prefix;N.result.witness=pair.second->witness;compare_dense_profile(*pair.first,N);}}
 certify_witness(p,coarse,RC,true);certify_witness(p,fine,RF,true);
 return {p,coarse,fine,t,K,cap,T.B,T.D,RC,RF,X,F.pilot.words,runs,positions,outside,ext,std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()};
}

#ifndef VDW_CUT_NO_MAIN
int main(int argc,char**argv){try{
 require(argc==8,"usage: cut p coarse fine pilot mode batch t");
 u32 p=parse(argv[1]),c=parse(argv[2]),f=parse(argv[3]),pd=parse(argv[4]),mode=parse(argv[5]),batch=parse(argv[6]),t=parse(argv[7]);
 require(mode<=4,"mode in 0..4");vdw_cut::Counts cut;
 auto O=shared_screened(p,c,f,t,pd,3,batch,false,mode,&cut);
 std::ostringstream tmp;auto* old=std::cout.rdbuf(tmp.rdbuf());print_shared(O);std::cout.rdbuf(old);
 std::string s=tmp.str();auto end=s.rfind('}');s.erase(end);
 std::cout<<s<<",\"cut_mode\":"<<mode<<",\"cut_words\":"<<cut.words<<",\"cut_queries\":"<<cut.queries<<",\"cut_rounds\":"<<cut.rounds<<",\"cut_input\":"<<cut.input<<",\"cut_output\":"<<cut.output<<",\"gpu_executed\":false}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
