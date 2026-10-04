#define VDW_ANCHOR_NO_MAIN
#include "../baseline/src/anchor_reference.cpp"
#include "periodic_bands.hpp"
#include "clipped_cover.hpp"
#include <tuple>
#include <array>
using Signature=std::tuple<u32,u32,u32,u32,u32,u32>;
struct ScheduleStats {
 u64 useful=0,old_slots=0,new_slots=0,row_tasks=0,row_bands=0;
 u64 old_rounds=0,new_rounds=0,old_active_warps=0,new_active_warps=0;
 u64 base_bank_rounds=0,base_conflicting_warps=0,duplicate_base_addresses=0;
 u32 max_bank=0;u64 periods[33]{};u64 swizzled_rounds[2]{};u32 swizzled_max[2]{};
};
// Simplified SM shared-memory bank model: 32 banks, distinct u32 addresses.
// Broadcasts to one address count once. This is NOT a measured instruction cost.
static u32 base_banks(const std::vector<vdw_periodic::Task>& tasks,u64&duplicates,u32 layout=0){
 std::array<std::vector<u32>,32> bank;
 for(auto x:tasks)if(x.valid){u32 a=x.word+1;if(layout>=1)a^=((x.word+1)>>5)&31;if(layout>=2)a^=((x.word+1)>>10)&31;bank[a&31].push_back(a);}
 u32 worst=0;for(auto&v:bank){std::sort(v.begin(),v.end());auto end=std::unique(v.begin(),v.end());u32 unique=end-v.begin();duplicates+=v.size()-unique;worst=std::max(worst,unique);}return worst;
}
static void extend_table(Table&T,vdw_clipped::Cover cv){
 require(vdw_clipped::certified(T.p,cv)&&cv.Q==T.L&&cv.B>=T.B,"extended cover domain");
 u32 old=T.B;T.B=cv.B;T.D=cv.D;T.c.resize(T.B+1);std::vector<u32>spf(T.B+1);
 for(u32 q=2;q<=T.B;++q)if(!spf[q])for(u32 x=q;x<=T.B;x+=q)if(!spf[x])spf[x]=q;
 std::vector<u64> roots(T.r);roots[0]=1;u64 z=power(T.g,(T.p-1)/T.r,T.p);for(u32 j=1;j<T.r;++j)roots[j]=roots[j-1]*z%T.p;
 for(u32 x=old+1;x<=T.B;++x){u32 q=spf[x];if(q<x)T.c[x]=(T.c[q]+T.c[x/q])%T.r;else{auto y=power(x,(T.p-1)/T.r,T.p);auto at=std::find(roots.begin(),roots.end(),y);require(at!=roots.end(),"extended character");T.c[x]=at-roots.begin();}}
 T.n=2*T.B+1;u32 np=0;while((1u<<np)<T.r)++np;T.planes.assign(np,std::vector<u32>((T.n+31)/32+2));
 for(u32 i=0;i<T.n;++i)if(i!=T.B){int cc=T.colour((i64)i-T.B);for(u32 j=0;j<np;++j)if((cc>>j)&1)T.planes[j][i>>5]|=1u<<(i&31);}
}
struct PeriodicExperiment {Experiment exp;ScheduleStats sched;u32 divisor=1;};
static PeriodicExperiment periodic_experiment(u32 p,u32 r,u32 pilot_d,bool audit=true,u32 divisor=1){
 require(prime(p)&&p>=7&&r>=2&&r<=9&&(p-1)%r==0&&divisor>=1&&divisor<=256,"experiment domain");
 auto begin=std::chrono::steady_clock::now();u32 t=std::max(3u,select_t(p,r).t);
 Table PT(p,r,t);AnchorStats initial;seed_table(PT,initial);AnchorStats pil=initial;
 if(pilot_d)pil=scan_bands(PT,t,t,std::min(pilot_d,PT.D),false,false,true,initial);
 u32 q=pil.maximum;require(q>=t,"pilot too short: dense fallback required");
 auto pc=power(pil.witness,(p-1)/r,p);for(u32 j=0;j<q;++j)require(power(pil.witness+j,(p-1)/r,p)==pc,"pilot witness");
 Table T(p,r,q);if(divisor>1)extend_table(T,vdw_clipped::from_denominator(p,q,std::max(1u,T.D/divisor)));Positive P(T);Wheel W2(T,2),W3(T,3),W6(T,6);auto inv=inverses(p,T.D);
 AnchorStats R;R.maximum=pil.maximum;R.witness=pil.witness;seed_table(T,R);
 ScheduleStats S;u32 spacing=q-t+1;std::vector<Signature> expected,emitted;
 auto start=std::chrono::steady_clock::now();
 for(u32 d=1;d<=T.D;++d){u32 g=std::gcd(d,6u);auto phase=vdw_anchor::choose(T.B,d,t,spacing,true,true);auto bands=vdw_anchor::bands(T.B,d,t,spacing,phase.phase);R.bands+=bands.size();
  auto traverse=[&](auto view,u32 c){
   auto pl=vdw_periodic::prepare(vdw_periodic::make_shapes(T.B,d,g,c,t,spacing,phase.phase));++S.row_tasks;
   for(auto side:{pl.a,pl.b})if(side.shape.middle)++S.periods[side.P];
   // Enumerated reference, also measures lane slots and stage envelopes.
   if(audit)for(auto b:bands){if(b.last<c)continue;u32 first=b.first>c?(b.first-c+g-1)/g:0,last=(b.last-c)/g;if(first>last)continue;++S.row_bands;
    u32 nw=(last>>5)-(first>>5)+1;S.old_slots+=32ull*((nw+31)/32);
    for(u32 j=0;j<nw;j+=32){u32 rounds=0;
     for(u32 lane=0;lane<32&&j+lane<nw;++lane){u32 word=(first>>5)+j+lane;vdw_next::Counters st;auto bits=vdw_next::filter_range(view,d/g,t,word,first,last,&st);rounds=std::max(rounds,(u32)st.stages);if(bits)expected.emplace_back(d,c,first,last,word,bits);}
     S.old_rounds+=rounds;++S.old_active_warps;
    }
   }
   S.new_slots+=32ull*((pl.total+31)/32);
   for(u32 rank=0;rank<pl.total;rank+=32){u32 rounds=0;std::vector<vdw_periodic::Task> tasks;
    for(u32 lane=0;lane<32&&rank+lane<pl.total;++lane){auto x=vdw_periodic::select(pl,rank+lane);require(x.valid,"packed selector");tasks.push_back(x);++S.useful;++R.words;
     vdw_next::Counters st;u32 bits=vdw_next::filter_range(view,d/g,t,x.word,x.first,x.last,&st);rounds=std::max(rounds,(u32)st.stages);R.stages+=st.stages;R.direct+=st.direct;R.signed_reads+=st.signed_reads;R.raw+=__builtin_popcount(bits);
     if(bits)emitted.emplace_back(d,c,x.first,x.last,x.word,bits);
     while(bits){u32 bit=__builtin_ctz(bits);bits&=bits-1;complete_anchor(T,t,g*(32*x.word+bit)+c,d,inv[d],R);}
    }
    S.new_rounds+=rounds;++S.new_active_warps;u32 banks=base_banks(tasks,S.duplicate_base_addresses);S.base_bank_rounds+=banks;S.max_bank=std::max(S.max_bank,banks);S.base_conflicting_warps+=banks>1;for(u32 a=1;a<=2;++a){u64 dup=0;u32 b=base_banks(tasks,dup,a);S.swizzled_rounds[a-1]+=b;S.swizzled_max[a-1]=std::max(S.swizzled_max[a-1],b);}
   }
  };
  if(g==1)traverse(P.view(),0);else{const Wheel&W=g==2?W2:g==3?W3:W6;for(u32 row=0;row<W.residues.size();++row)traverse(W.view(row),W.residues[row]);}
 }
 if(audit){std::sort(expected.begin(),expected.end());std::sort(emitted.begin(),emitted.end());require(expected==emitted,"candidate masks differ from previous enumerator");}
 double mainseconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();bool checked=p<=2000000;
 if(checked){NewResult N;N.result.M=R.maximum;N.result.prefix=R.prefix;N.result.witness=R.witness;compare_dense_profile(T,N);}
 auto ch=power(R.witness,(p-1)/r,p);for(u32 j=0;j<R.maximum;++j)require(power(R.witness+j,(p-1)/r,p)==ch,"direct witness");
 if(R.witness>1)require(power(R.witness-1,(p-1)/r,p)!=ch,"left neighbour");if((u64)R.witness+R.maximum<p)require(power(R.witness+R.maximum,(p-1)/r,p)!=ch,"right neighbour");
 return {{p,r,t,std::min(pilot_d,PT.D),q,T.B,T.D,pil,R,checked,std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count(),mainseconds},S,divisor};
}
static void print_periodic(const PeriodicExperiment&E){auto&S=E.sched;auto&P=E.exp;auto&R=P.main;
 std::cout<<std::setprecision(10)<<"{\"status\":\"PASS\",\"p\":"<<P.p<<",\"denominator_divisor\":"<<E.divisor<<",\"r\":"<<P.r<<",\"t\":"<<P.t<<",\"Q\":"<<P.q<<",\"B\":"<<P.B<<",\"D\":"<<P.D<<",\"maximum\":"<<R.maximum<<",\"prefix\":"<<R.prefix<<",\"witness\":"<<R.witness<<",\"pilot_words\":"<<P.pilot.words
 <<",\"main_words\":"<<R.words<<",\"main_stages\":"<<R.stages<<",\"raw_bits\":"<<R.raw<<",\"primitive_anchors\":"<<R.primitives<<",\"extensions\":"<<R.extensions<<",\"oracle_queries\":"<<R.queries
 <<",\"row_step_tasks\":"<<S.row_tasks<<",\"row_bands\":"<<S.row_bands<<",\"old_lane_slots\":"<<S.old_slots<<",\"packed_lane_slots\":"<<S.new_slots<<",\"old_stage_envelope\":"<<S.old_rounds<<",\"packed_stage_envelope\":"<<S.new_rounds
 <<",\"old_active_warps\":"<<S.old_active_warps<<",\"packed_active_warps\":"<<S.new_active_warps<<",\"base_bank_rounds_model\":"<<S.base_bank_rounds<<",\"base_conflicting_warps_model\":"<<S.base_conflicting_warps<<",\"max_base_bank_multiplicity_model\":"<<S.max_bank<<",\"duplicate_base_addresses\":"<<S.duplicate_base_addresses<<",\"swizzled_base_bank_rounds_model\":["<<S.swizzled_rounds[0]<<","<<S.swizzled_rounds[1]<<"],\"swizzled_max_bank_model\":["<<S.swizzled_max[0]<<","<<S.swizzled_max[1]<<"],\"period_counts\":[";
 for(u32 k=1;k<=32;k*=2)std::cout<<(k>1?",":"")<<S.periods[k];
 std::cout<<"],\"candidate_masks_equal\":true,\"dense_checked\":"<<(P.dense_checked?"true":"false")<<",\"GPU_executed\":false,\"cpu_seconds_with_two_filter_paths_and_models\":"<<P.seconds<<"}\n";
}
#ifndef VDW_PERIODIC_NO_MAIN
int main(int argc,char**argv){try{require((argc==5||argc==6)&&std::string(argv[1])=="--profile","usage: --profile PRIME R PILOT_D [D_DIVISOR]");print_periodic(periodic_experiment(parse(argv[2]),parse(argv[3]),parse(argv[4]),true,argc==6?parse(argv[5]):1));}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
