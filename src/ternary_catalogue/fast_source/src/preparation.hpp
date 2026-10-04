#pragma once
#include "mont32.hpp"
#define VDW_TARGET_NO_MAIN
#include "../vendor/research/src/target_reference.cpp"
#include "completion_core.hpp"
#include <memory>
#include <iomanip>
namespace kit {
struct HostPlanes {std::vector<u32> data;u32 stride=0,np=0,flip=0;};
struct Prep {
 std::unique_ptr<Table> fine,coarse;
 HostPlanes cp[4],fp[4];std::vector<u32> inv;
 AnchorStats rc,rf;u32 t=0,cap=0,D=0,B=0,selected=0,p=0,c=0,f=0;
 double seconds=0;
 Prep(u32 pp,u32 cc,u32 ff,u32 pilot):p(pp),c(cc),f(ff){
  auto start=std::chrono::steady_clock::now();
  require(prime(p)&&c>=2&&f>=c&&f<=9&&f%c==0&&(p-1)%f==0,"group domain");
  Prepared seed(p,f,pilot);u32 K=seed.pilot.maximum;cap=K+3;t=std::max(3u,(K+1)/2);
  require(K>=t&&cap<p&&cap<=256,"witnessed seed too small / unsupported cover");
  fine=std::make_unique<Table>(p,f,cap);
  extend_table(*fine,vdw_clipped::from_denominator(p,cap,std::max(1u,fine->D/2)));
  B=fine->B;D=fine->D;coarse=std::make_unique<Table>(*fine);
  coarse->r=c;coarse->sigma%=c;for(auto &x:coarse->c)x%=c;
  u32 np=0;while((1u<<np)<c)++np;coarse->planes.assign(np,{});
  rf.maximum=rc.maximum=K;rf.witness=rc.witness=seed.pilot.witness;
  seed_table(*coarse,rc);seed_table(*fine,rf);
  certify_witness(p,c,rc,false);certify_witness(p,f,rf,false);
  require(rc.maximum>=rf.maximum,"seed refinement invariant");
  auto fill=[](Table &T,HostPlanes*o){Positive pos(T);o[0]={pos.data,pos.stride,pos.np,pos.flip};
   int ix=1;for(u32 g:{2u,3u,6u}){Wheel w(T,g);o[ix++]={w.data,w.stride,w.np,w.flip};}};
  fill(*coarse,cp);fill(*fine,fp);inv=inverses(p,D);
#if VDW_COMPLETION_MONT
  {u32 one=(1ull<<32)%p;for(auto &x:inv)x=(u64)x*one%p;}
#endif
  selected=c<f?minimal_planes(c,f):0;
  seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 }
 State seed_state()const{State s{};s.best_c=vdw_next::pack_best(rc.maximum,rc.witness);s.best_f=vdw_next::pack_best(rf.maximum,rf.witness);s.cap=cap;s.D=D;s.t=t;s.coarse=c;s.fine=f;return s;}
 vdw_next::CompletionView completion()const{vdw_next::CompletionView v{p,f,B,t,D,fine->c.data(),inv.data()};
#if VDW_COMPLETION_MONT
  v.ni=kit::ninverse(p);
#endif
  v.sigma_value=((p-1)/2)%f;return v;}
};
inline std::vector<std::pair<u32,u32>> groups(u32 p,const std::string &scope){
 std::vector<std::pair<u32,u32>> v;
 if(scope=="24"){require((p-1)%4==0,"24 needs p=1 mod4");return {{2,4}};}
 if(scope=="36"){require((p-1)%6==0,"36 needs p=1 mod6");return {{3,6}};}
 require(scope=="core"||scope=="all","scope must be core, all, 24 or 36");
 if((p-1)%4==0)v.push_back({2,4});else v.push_back({2,2});
 if((p-1)%3==0)v.push_back({3,6});
 if(scope=="all")for(u32 r:{5u,7u,8u,9u})if((p-1)%r==0)v.push_back({r,r});
 return v;
}
inline u32 requested_mask(u32 p,const std::string&s){u32 m=0;for(auto [c,f]:groups(p,s)){m|=1u<<(c-2);m|=1u<<(f-2);}return m;}
inline void dense_host(u32 p,u32 mask,u32 *m,u32 *pre,u32*witness){
 require(p<=2000000&&prime(p),"dense CPU oracle domain");u32 g=root(p);std::vector<unsigned short> idx(p);u64 x=1;
 for(u32 i=0;i<p-1;++i){idx[x]=i%2520;x=x*g%p;}
 for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2))){u32 run=0,previous=UINT32_MAX;m[r]=0;pre[r]=1;
  for(u32 j=1;j<p;++j){u32 col=idx[j]%r;run=col==previous?run+1:1;previous=col;if(run>m[r]){m[r]=run;witness[r]=j-run+1;}}
  while(pre[r]+1<p&&idx[pre[r]+1]%r==0)++pre[r];
 }
}
inline State host_integrated(Prep &p,u32 epoch,bool target){
 State s=p.seed_state();auto cv=p.completion();
 Positive C(*p.coarse),F(*p.fine);Wheel C2(*p.coarse,2),C3(*p.coarse,3),C6(*p.coarse,6),F2(*p.fine,2),F3(*p.fine,3),F6(*p.fine,6);
 for(u32 lo=1;lo<=p.D;lo+=epoch){take_snapshot(s,lo,std::min(p.D+1,lo+epoch));
  for(u32 d=lo;d<s.hi;++d){u32 g=std::gcd(d,6u),space=s.q-p.t+1,phase=((p.t-1)/2)%space,tag=g==1?0:g==2?1:g==3?2:3;
   auto row=[&](auto c,auto f,u32 residue){auto plan=vdw_periodic::prepare(vdw_periodic::make_shapes(p.B,d,g,residue,p.t,space,phase));
    for(u32 rank=0;rank<plan.total;++rank){auto task=vdw_periodic::select(plan,rank);u32 bits=vdw_next::filter_range(c,d/g,p.t,task.word,task.first,task.last);
     vdw_target::Masks mm;
     if(target&&p.c<p.f)mm=vdw_target::split(c,f,p.selected,d/g,p.t,s.qc,s.qf,task.word,bits);
     else mm={vdw_cut::screen(c,d/g,p.t,s.q,task.word,bits,1,nullptr),0};
     if(mm.coarse|mm.fine){u64 bc=s.snap_c,bf=s.snap_f;require(kit::complete(cv,p.c,{d|(tag<<30),g*32*task.word+residue,mm.coarse,mm.fine},bc,bf),"CPU integrated completion");s.best_c=std::max<unsigned long long>(s.best_c,bc);s.best_f=std::max<unsigned long long>(s.best_f,bf);++s.count;}
    }
   };
   if(g==1)row(C.view(),F.view(),0);else{auto&cw=g==2?C2:g==3?C3:C6;auto&fw=g==2?F2:g==3?F3:F6;for(u32 i=0;i<cw.residues.size();++i)row(cw.view(i),fw.view(i),cw.residues[i]);}
  }++s.epochs_done;
 }return s;
}
} // namespace kit
