#ifndef KIT_REUSE_EVENTS
#define KIT_REUSE_EVENTS 0
#endif
#include "gpu_pipeline.cuh"
#include "small_table.cuh"
#if KIT_EXACT_TABLE_SEED
#include "seed_boundary.cuh"
#endif
#include <cuda_profiler_api.h>
#include <sstream>
#include <cstdio>
#ifndef KIT_CHARACTER_ROOT
#define KIT_CHARACTER_ROOT 0
#endif
// Only a primitive m-th root is required. Changing that root multiplies all
// colour indices by a unit modulo m, preserving equality, divisor refinements,
// zero prefixes, reflection and the canonical smallest witnesses.
static u32 character_root(u32 p,u32 m){
 require(m>=2&&m<=2520&&(p-1)%m==0,"character root domain");
 if(m==2)return p-1;
 u32 n=m;std::vector<u32>factors;for(u32 q=2;q*q<=n;++q)if(n%q==0){factors.push_back(q);do n/=q;while(n%q==0);}if(n>1)factors.push_back(n);
 for(u32 a=2;a<p;++a){u32 z=(u32)power(a,(p-1)/m,p);bool ok=z!=0;for(u32 q:factors)if(power(z,m/q,p)==1){ok=false;break;}if(ok)return z;}
 throw std::runtime_error("primitive character root not found");
}
extern "C" void baseline_init(unsigned);
extern "C" int baseline_profile(unsigned,unsigned,unsigned*,unsigned*,int);
static void ck(cudaError_t e){if(e!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(e));}
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a){return std::chrono::duration<double,std::milli>(Clock::now()-a).count();}
struct Options {std::string engine="candidate",scope="all",memory="shared",prep="gpu";unsigned grid_policy=0,cover_q=16,denominator_divisor=2,generator_blocks=2,completion_blocks=2,single3_q=0,single6_q=0,other_q=0,group36_q=0,single4_q=0,pair_mode=0,gpu_pilot=0,dynamic_t=0,t_mode=0,pilot_threads=0,shared_min_resident=0,threads=256,epoch=1024,pilot=64,repeats=1,capacity=262144;bool graph=false,target=true,profile=false;std::vector<u32> primes;};
static Options options(int argc,char**argv){Options o;
 for(int i=1;i<argc;++i){std::string k=argv[i];if(k=="--profile"){o.profile=true;continue;}if(i+1==argc)throw std::runtime_error("missing option value");std::string v=argv[++i];
  if(k=="--single3-q")o.single3_q=parse(v.c_str());else if(k=="--single6-q")o.single6_q=parse(v.c_str());else if(k=="--other-q")o.other_q=parse(v.c_str());else if(k=="--group36-q")o.group36_q=parse(v.c_str());else if(k=="--single4-q")o.single4_q=parse(v.c_str());else if(k=="--pair-mode")o.pair_mode=parse(v.c_str());else if(k=="--t-mode")o.t_mode=parse(v.c_str());else if(k=="--pilot-threads")o.pilot_threads=parse(v.c_str());else if(k=="--gpu-pilot")o.gpu_pilot=parse(v.c_str());else if(k=="--dynamic-t")o.dynamic_t=parse(v.c_str());else if(k=="--grid-policy")o.grid_policy=parse(v.c_str());else if(k=="--cover-q")o.cover_q=parse(v.c_str());else if(k=="--denominator-divisor")o.denominator_divisor=parse(v.c_str());else if(k=="--generator-blocks")o.generator_blocks=parse(v.c_str());else if(k=="--completion-blocks")o.completion_blocks=parse(v.c_str());else if(k=="--shared-min-resident")o.shared_min_resident=parse(v.c_str());else if(k=="--engine")o.engine=v;else if(k=="--scope")o.scope=v;else if(k=="--memory")o.memory=v;else if(k=="--prep")o.prep=v;
  else if(k=="--threads")o.threads=parse(v.c_str());else if(k=="--epoch")o.epoch=parse(v.c_str());else if(k=="--pilot")o.pilot=parse(v.c_str());
  else if(k=="--repeat")o.repeats=parse(v.c_str());else if(k=="--capacity")o.capacity=parse(v.c_str());
  else if(k=="--graph"){require(v=="0"||v=="1","graph boolean");o.graph=v=="1";}
  else if(k=="--target"){require(v=="0"||v=="1","target boolean");o.target=v=="1";}
  else if(k=="--primes"){std::istringstream s(v);std::string q;while(std::getline(s,q,','))o.primes.push_back(parse(q.c_str()));}
  else throw std::runtime_error("unknown option: "+k);
 }
 require(!o.primes.empty()&&o.primes.size()<=128,"1..128 primes required");
 require(o.engine=="candidate"||o.engine=="baseline"||o.engine=="cpu","engine");
 require(o.prep=="cpu"||o.prep=="gpu","prep");require(o.memory=="shared"||o.memory=="l2","memory");require(o.threads>=128&&o.threads<=1024&&o.threads%32==0,"threads");
 require(o.grid_policy<=3,"grid policy");
 require(o.gpu_pilot<=8192&&o.dynamic_t<=1,"bounded GPU pilot / dynamic threshold");
 require(o.pair_mode<=3,"pair splitting mode");
 require(!o.single4_q||(o.single4_q>=3&&o.single4_q<=o.cover_q&&(o.pair_mode&1)&&o.prep=="gpu"),"single4 cover requires GPU split24 and 3<=Q4<=Q");
 require(!o.group36_q||(o.group36_q>=3&&o.group36_q<=o.cover_q&&!(o.pair_mode&2)&&o.prep=="gpu"),"36 cover requires coupled GPU group and 3<=Q36<=Q");
 require(!o.other_q||(o.other_q>=3&&o.other_q<=o.cover_q&&o.prep=="gpu"),"singleton cover requires GPU prep and 3<=Qr<=Q");
 require((!o.single3_q||(o.single3_q>=3&&o.single3_q<=o.cover_q))&&(!o.single6_q||(o.single6_q>=3&&o.single6_q<=o.cover_q)),"split36 cover range");
 require((!o.single3_q&&!o.single6_q)||((o.pair_mode&2)&&o.prep=="gpu"),"split36 covers require GPU split36");
 require(o.t_mode<=13&&(o.pilot_threads==0||o.pilot_threads==128||o.pilot_threads==256||o.pilot_threads==512||o.pilot_threads==1024),"anchor mode / pilot threads");
 require(o.cover_q>=3&&o.cover_q<=256&&o.denominator_divisor>=1&&o.denominator_divisor<=16,"cover options");
 require(o.prep=="gpu"||(o.cover_q==16&&o.denominator_divisor==2),"cover tuning is GPU-only");
 require(o.generator_blocks>=1&&o.generator_blocks<=8&&o.completion_blocks>=1&&o.completion_blocks<=64&&o.shared_min_resident<=3,"bounded experiments");
 require(o.epoch>=32&&o.epoch<=65536&&o.pilot<=512&&o.repeats>=1&&o.repeats<=16&&o.capacity>=1&&o.capacity<=4194304,"bounded options");
 for(u32 p:o.primes){require(p>=1000&&prime(p),"prime domain [1000,2^32)");kit::requested_mask(p,o.scope);}return o;
}
// Disjoint denominator intervals retain complete coverage exactly once.
static std::vector<std::pair<u32,u32>> gpu_epochs(u32 D,u32 epoch,u32 pilot){
 std::vector<std::pair<u32,u32>> out;u32 lo=1;
 if(pilot&&pilot<D){out.push_back({1,pilot+1});lo=pilot+1;}
 for(;lo<=D;lo+=epoch)out.push_back({lo,std::min(D+1,lo+epoch)});
 return out;
}
struct Workspace {
 cudaStream_t stream{};cudaEvent_t timer_begin{},timer_end{};kit::State *state{},*hs{};kit::Candidate*out{};unsigned char*colour{},*hcolour{};u32 *inverse{},*hinverse{};
 u32 *cp[4]{},*fp[4]{},*hcp[4]{},*hfp[4]{};size_t plan_words=150000;unsigned capacity;cudaDeviceProp prop{};unsigned short*base{};unsigned*base_error{};u32*roots{},*labels{},*hroots{},*hlabels{};
 Workspace(unsigned cap):capacity(cap){ck(cudaGetDeviceProperties(&prop,0));require(prop.major==12&&prop.minor==0,"this kit is qualified only for SM120");
  ck(cudaMalloc(&base,1048577*2ull));ck(cudaMalloc(&base_error,4));ck(cudaMalloc(&roots,2560*4));ck(cudaMalloc(&labels,2560*4));ck(cudaMallocHost(&hroots,2560*4));ck(cudaMallocHost(&hlabels,2560*4));
  ck(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));if constexpr(KIT_REUSE_EVENTS){ck(cudaEventCreate(&timer_begin));ck(cudaEventCreate(&timer_end));}ck(cudaMalloc(&state,sizeof(*state)));ck(cudaMallocHost(&hs,sizeof(*hs)));
  ck(cudaMalloc(&out,(size_t)cap*sizeof(*out)));ck(cudaMalloc(&colour,1048577));ck(cudaMallocHost(&hcolour,1048577));
  ck(cudaMalloc(&inverse,1048577*4ull));ck(cudaMallocHost(&hinverse,1048577*4ull));
  for(int i=0;i<4;++i){ck(cudaMalloc(&cp[i],plan_words*4));ck(cudaMalloc(&fp[i],plan_words*4));ck(cudaMallocHost(&hcp[i],plan_words*4));ck(cudaMallocHost(&hfp[i],plan_words*4));}
 }
 ~Workspace(){cudaStreamSynchronize(stream);cudaFree(base);cudaFree(base_error);cudaFree(roots);cudaFree(labels);cudaFreeHost(hroots);cudaFreeHost(hlabels);cudaStreamSynchronize(stream);for(int i=0;i<4;++i){cudaFree(cp[i]);cudaFree(fp[i]);cudaFreeHost(hcp[i]);cudaFreeHost(hfp[i]);}cudaFree(state);cudaFree(out);cudaFree(colour);cudaFree(inverse);cudaFreeHost(hs);cudaFreeHost(hcolour);cudaFreeHost(hinverse);if constexpr(KIT_REUSE_EVENTS){cudaEventDestroy(timer_begin);cudaEventDestroy(timer_end);}cudaStreamDestroy(stream);}
};
static int grid_for_epoch(const Options&o,const Workspace&w,u32 D,u32 lo,u32 hi){
 u32 maximum=w.prop.multiProcessorCount*o.generator_blocks;
 if(!o.grid_policy)return maximum;
 u32 cycle0=(lo-1)/6,cycles=(hi-1+5)/6-cycle0,n=cycles*2,slots;
#if KIT_ADAPTIVE_CHUNKS
 u32 j0=cycle0<KIT_LOW_CYCLES?std::min(n,(KIT_LOW_CYCLES-cycle0)*2):0,half=(D+6*KIT_MID_DEN_DIV-1)/(6*KIT_MID_DEN_DIV);
 u32 j1=std::max(j0,cycle0<half?std::min(n,(half-cycle0)*2):0);
 slots=(j0<<KIT_LOW_CHUNK_SHIFT)+((j1-j0)<<KIT_MID_CHUNK_SHIFT)+((n-j1)<<KIT_HIGH_CHUNK_SHIFT);
#else
 slots=n<<KIT_CHUNK_SHIFT;
#endif
 u32 warps=o.threads/32,blocks=(slots+warps-1)/warps;
 u32 floor=o.grid_policy==1?w.prop.multiProcessorCount:o.grid_policy==2?(w.prop.multiProcessorCount+1)/2:1;
 return std::min(maximum,std::max(floor,blocks));
}
struct Metrics{std::string group_diagnostics;double prep=0,transfer=0,gpu_prepare=0,setup=0,queue=0,verify=0;unsigned groups=0,fallbacks=0;std::string resources;};
static void run_group(Workspace&w,kit::Prep &p,const Options&o,u32*m,u32*pre,u32*wit,Metrics&stats){
 auto setupstart=Clock::now();kit::GInput c[4],f[4];bool use[4];kit::Resources res[4];u32 gs[4]={1,2,3,6};
 for(int i=0;i<4;++i){require(p.cp[i].data.size()<=w.plan_words&&p.fp[i].data.size()<=w.plan_words,"plane pool capacity");
  std::copy(p.cp[i].data.begin(),p.cp[i].data.end(),w.hcp[i]);std::copy(p.fp[i].data.begin(),p.fp[i].data.end(),w.hfp[i]);
  c[i]={w.cp[i],p.cp[i].stride,p.cp[i].np,p.B,p.cp[i].flip};f[i]={w.fp[i],p.fp[i].stride,p.fp[i].np,p.B,p.fp[i].flip};
  use[i]=o.memory=="shared"&&p.cp[i].data.size()*4<=(w.prop.sharedMemPerBlockOptin-kit::PLAN_CACHE_BYTES);
  ck(kit::launch_group(gs[i],c[i],f[i],w.state,p.selected,o.target,w.out,w.capacity,w.prop.multiProcessorCount*o.generator_blocks,o.threads,use[i],w.stream,&res[i]));
  std::ostringstream r;r<<"{\"r\":"<<p.c<<",\"fine\":"<<p.f<<",\"wheel\":"<<gs[i]<<",\"registers\":"<<res[i].registers<<",\"local_bytes\":"<<res[i].local_bytes<<",\"resident_blocks\":"<<res[i].active_blocks<<",\"shared_bytes\":"<<res[i].shared_bytes<<",\"shared\":"<<(use[i]?"true":"false")<<"}";
  if(!stats.resources.empty())stats.resources+=",";stats.resources+=r.str();
 }
 *w.hs=p.seed_state();std::copy(p.fine->c.begin(),p.fine->c.end(),w.hcolour);std::copy(p.inv.begin(),p.inv.end(),w.hinverse);
 stats.setup+=ms(setupstart);cudaEvent_t a{},b{};if constexpr(KIT_REUSE_EVENTS){a=w.timer_begin;b=w.timer_end;}else{ck(cudaEventCreate(&a));ck(cudaEventCreate(&b));}
 ck(cudaEventRecord(a,w.stream));
 for(int i=0;i<4;++i){ck(cudaMemcpyAsync(w.cp[i],w.hcp[i],p.cp[i].data.size()*4,cudaMemcpyHostToDevice,w.stream));ck(cudaMemcpyAsync(w.fp[i],w.hfp[i],p.fp[i].data.size()*4,cudaMemcpyHostToDevice,w.stream));}
 ck(cudaMemcpyAsync(w.colour,w.hcolour,p.fine->c.size(),cudaMemcpyHostToDevice,w.stream));ck(cudaMemcpyAsync(w.inverse,w.hinverse,p.inv.size()*4,cudaMemcpyHostToDevice,w.stream));
 ck(cudaMemcpyAsync(w.state,w.hs,sizeof(*w.hs),cudaMemcpyHostToDevice,w.stream));ck(cudaEventRecord(b,w.stream));ck(cudaEventSynchronize(b));float x;ck(cudaEventElapsedTime(&x,a,b));stats.transfer+=x;
 vdw_next::CompletionView cv{p.p,p.f,p.B,p.t,p.D,w.colour,w.inverse};
#if VDW_COMPLETION_MONT
 cv.ni=kit::ninverse(p.p);
#endif
 cv.sigma_value=((p.p-1)/2)%p.f;

 auto chain=[&](){for(u32 lo=1;lo<=p.D;lo+=o.epoch){
   kit::epoch_snapshot<<<1,1,0,w.stream>>>(w.state,lo,std::min(p.D+1,lo+o.epoch));
   for(int i=0;i<4;++i)ck(kit::launch_group(gs[i],c[i],f[i],w.state,p.selected,o.target,w.out,w.capacity,grid_for_epoch(o,w,p.D,lo,std::min(p.D+1,lo+o.epoch)),o.threads,use[i],w.stream));
   kit::complete_candidates<<<w.prop.multiProcessorCount*o.completion_blocks,128,0,w.stream>>>(cv,w.state,w.out,w.capacity);
   kit::end_epoch<<<1,1,0,w.stream>>>(w.state,w.capacity);ck(cudaPeekAtLastError());
  }};
 cudaGraph_t graph{};cudaGraphExec_t exec{};
 auto gs_start=Clock::now();
 if(o.graph){ck(cudaStreamBeginCapture(w.stream,cudaStreamCaptureModeThreadLocal));chain();ck(cudaStreamEndCapture(w.stream,&graph));ck(cudaGraphInstantiate(&exec,graph,0));}
 stats.setup+=ms(gs_start);
 ck(cudaEventRecord(a,w.stream));if(o.graph)ck(cudaGraphLaunch(exec,w.stream));else chain();
 ck(cudaMemcpyAsync(w.hs,w.state,sizeof(*w.hs),cudaMemcpyDeviceToHost,w.stream));ck(cudaEventRecord(b,w.stream));ck(cudaEventSynchronize(b));ck(cudaEventElapsedTime(&x,a,b));stats.queue+=x;
 require(!w.hs->failed,"GPU error/overflow: result refused");require(w.hs->epochs_done==(p.D+o.epoch-1)/o.epoch,"incomplete denominator coverage");
 m[p.c]=w.hs->best_c>>32;m[p.f]=w.hs->best_f>>32;wit[p.c]=0xffffffffu-(u32)w.hs->best_c;wit[p.f]=0xffffffffu-(u32)w.hs->best_f;pre[p.c]=p.rc.prefix;pre[p.f]=p.rf.prefix;
 if(exec)ck(cudaGraphExecDestroy(exec));if(graph)ck(cudaGraphDestroy(graph));if constexpr(!KIT_REUSE_EVENTS){ck(cudaEventDestroy(a));ck(cudaEventDestroy(b));}++stats.groups;
}
static vdw_clipped::Cover prime_cover(u32 p,u32 Q,const Options&o){
 u64 area=(u64)(Q-1)*p,B=(u64)std::sqrt((long double)area);
 while((B+1)*(B+1)<=area)++B;while(B*B>area)--B;++B;
 u32 D=std::max(1u,(u32)(2*B/(Q-1))/o.denominator_divisor);
 while((u64)(Q-1)*D*D>4ull*p)--D;
 auto cover=vdw_clipped::from_denominator(p,Q,D);
 return cover;
}
// Every group gets one independently certified cover. Fine singletons may
// use a smaller table; common character preparation still covers the largest B.
static vdw_clipped::Cover group_cover(u32 p,u32 c,u32 f,vdw_clipped::Cover base,const Options&o){
 u32 Q=0;
 if(c==f){if(f==3)Q=o.single3_q;if(f==4)Q=o.single4_q;
  if(f>=5)Q=o.other_q;if(f==6&&o.single6_q)Q=o.single6_q;
 }else if(c==3&&f==6)Q=o.group36_q;
 return Q?prime_cover(p,Q,o):base;
}
static vdw_clipped::Cover prepare_gpu_table(Workspace&w,u32 p,u32 mask,const Options&o,Metrics&s){
 auto a=Clock::now();auto cover=prime_cover(p,o.cover_q,o);
 u32 inverse_D=cover.D;
 for(auto [c,f]:kit::groups(p,o.scope)){
  bool split=c<f&&((c==2&&(o.pair_mode&1))||(c==3&&(o.pair_mode&2)));
  auto ca=group_cover(p,c,split?c:f,cover,o),cb=split?group_cover(p,f,f,cover,o):ca;
  require(ca.B<=cover.B&&cb.B<=cover.B,"group table beyond prepared base");
  inverse_D=std::max(inverse_D,std::max(ca.D,cb.D));
 }
 u32 m=1;for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2)))m=std::lcm(m,r);
 u32 ni=kit::ninverse(p),one=(u32)((1ull<<32)%p),R2=(u64)one*one%p;
#if KIT_CHARACTER_ROOT
 u32 z=character_root(p,m);
#else
 u32 z=(u32)power(root(p),(p-1)/m,p);
#endif
 u32 cur=1;std::vector<std::pair<u32,u32>> v(m);
 for(u32 j=0;j<m;++j){v[j]={kit::mont(cur,R2,p,ni),j};cur=(u64)cur*z%p;}std::sort(v.begin(),v.end());
 for(u32 j=0;j<m;++j){w.hroots[j]=v[j].first;w.hlabels[j]=v[j].second;}
 s.prep+=ms(a);cudaEvent_t b{},e{};if constexpr(KIT_REUSE_EVENTS){b=w.timer_begin;e=w.timer_end;}else{ck(cudaEventCreate(&b));ck(cudaEventCreate(&e));}ck(cudaEventRecord(b,w.stream));
 ck(cudaMemcpyAsync(w.roots,w.hroots,m*4,cudaMemcpyHostToDevice,w.stream));ck(cudaMemcpyAsync(w.labels,w.hlabels,m*4,cudaMemcpyHostToDevice,w.stream));ck(cudaMemsetAsync(w.base_error,0,4,w.stream));
 kit::small_table<<<(cover.B+256)/256,256,0,w.stream>>>(p,cover.B,m,ni,R2,one,w.roots,w.labels,w.base,w.base_error);
 kit::inverse_table<<<(inverse_D+256)/256,256,0,w.stream>>>(p,inverse_D,ni,R2,one,w.inverse);ck(cudaPeekAtLastError());
 ck(cudaEventRecord(e,w.stream));ck(cudaEventSynchronize(e));float dt;ck(cudaEventElapsedTime(&dt,b,e));s.gpu_prepare+=dt;
 if constexpr(!KIT_REUSE_EVENTS){ck(cudaEventDestroy(b));ck(cudaEventDestroy(e));}return cover;
}
static void gpu_prepared_group(Workspace&w,u32 p,u32 coarse,u32 fine,vdw_clipped::Cover cover,const Options&o,u32*m,u32*pre,u32*wit,Metrics&s){
 auto host_start=Clock::now();u32 np=0,nf=0;while((1u<<np)<coarse)++np;while((1u<<nf)<fine)++nf;
 auto flip=[&](u32 r){u32 two=1;while(r%(2*two)==0)two*=2;return ((p-1)/2)%r?two/2:0u;};
 kit::GInput c[4],f[4];bool shared[4];kit::Resources resources[4];u32 gs[4]={1,2,3,6};
 u32 selected=coarse<fine?minimal_planes(coarse,fine):0;
 for(int i=0;i<4;++i){u32 stride=(cover.B/gs[i]+32)/32+3,rows=i<2?1:2;
  c[i]={w.cp[i],stride,np,cover.B,flip(coarse)};f[i]={w.fp[i],stride,nf,cover.B,flip(fine)};
#if KIT_SINGLETON_FAST
  if(coarse==fine)f[i].data=w.cp[i];
#endif
  shared[i]=o.memory=="shared"&&(size_t)rows*np*stride*4<=(w.prop.sharedMemPerBlockOptin-kit::PLAN_CACHE_BYTES);
  ck(kit::launch_group(gs[i],c[i],f[i],w.state,selected,o.target,w.out,w.capacity,w.prop.multiProcessorCount*o.generator_blocks,o.threads,shared[i],w.stream,&resources[i]));
  if(shared[i]&&resources[i].active_blocks<(int)o.shared_min_resident){shared[i]=false;ck(kit::launch_group(gs[i],c[i],f[i],w.state,selected,o.target,w.out,w.capacity,w.prop.multiProcessorCount*o.generator_blocks,o.threads,shared[i],w.stream,&resources[i]));}
  std::ostringstream r;r<<"{\"r\":"<<coarse<<",\"fine\":"<<fine<<",\"wheel\":"<<gs[i]<<",\"registers\":"<<resources[i].registers<<",\"local_bytes\":"<<resources[i].local_bytes<<",\"resident_blocks\":"<<resources[i].active_blocks<<",\"shared_bytes\":"<<resources[i].shared_bytes<<",\"shared\":"<<(shared[i]?"true":"false")<<"}";
  if(!s.resources.empty())s.resources+=",";s.resources+=r.str();
 }
 *w.hs={};w.hs->best_c=w.hs->best_f=vdw_next::pack_best(1,1);w.hs->cap=cover.Q;w.hs->D=cover.D;w.hs->t=3;w.hs->coarse=coarse;w.hs->fine=fine;
 vdw_next::CompletionView cv{p,fine,cover.B,3,cover.D,w.colour,w.inverse};
#if VDW_COMPLETION_MONT
 cv.ni=kit::ninverse(p);
#endif
 cv.sigma_value=((p-1)/2)%fine;
auto epochs=gpu_epochs(cover.D,o.epoch,o.gpu_pilot);
s.setup+=ms(host_start);
 auto chain=[&](){
  ck(cudaMemcpyAsync(w.state,w.hs,sizeof(*w.hs),cudaMemcpyHostToDevice,w.stream));
  kit::fine_labels<<<(cover.B+256)/256,256,0,w.stream>>>(cover.B,fine,w.base,w.colour);
  for(int i=0;i<4;++i){u32 rows=i<2?1:2;
   #if KIT_WARP_PACK
   kit::pack_pair<<<(rows*c[i].stride+7)/8,256,0,w.stream>>>(cover.B,coarse,fine,gs[i],np,nf,c[i].stride,w.colour,w.cp[i],
#if KIT_SINGLETON_FAST
    coarse==fine?w.cp[i]:w.fp[i]
#else
    w.fp[i]
#endif
   );
#else
   kit::pack_planes<<<(rows*np*c[i].stride+255)/256,256,0,w.stream>>>(cover.B,coarse,gs[i],np,c[i].stride,w.colour,w.cp[i]);
   kit::pack_planes<<<(rows*nf*f[i].stride+255)/256,256,0,w.stream>>>(cover.B,fine,gs[i],nf,f[i].stride,w.colour,w.fp[i]);
#endif
  }
  u32 segments=(cover.B+1023)/1024;
#if KIT_SINGLETON_SEED
  if(coarse==fine)kit::seed_single<<<(segments+3)/4,128,0,w.stream>>>(p,cover.B,w.colour,w.state);else
#endif
  kit::seed_runs<<<KIT_WARP_SEED?(segments+3)/4:(segments+127)/128,128,0,w.stream>>>(p,cover.B,coarse,fine,w.colour,w.state);
#if KIT_EXACT_TABLE_SEED
  kit::seed_boundaries<<<((cover.B-1)/1024+128)/128,128,0,w.stream>>>(cv,w.state);
#endif
  kit::finish_seed<<<1,1,0,w.stream>>>(cv,w.state,w.base_error,o.t_mode);
  for(auto [lo,hi]:epochs){
   kit::epoch_snapshot<<<1,1,0,w.stream>>>(w.state,lo,hi,o.dynamic_t&&lo>1,o.t_mode);
   Options launch=o;if(lo==1&&o.gpu_pilot&&o.pilot_threads)launch.threads=o.pilot_threads;
   for(int i=0;i<4;++i)ck(kit::launch_group(gs[i],c[i],f[i],w.state,selected,o.target,w.out,w.capacity,grid_for_epoch(launch,w,cover.D,lo,hi),launch.threads,shared[i],w.stream));
   kit::complete_candidates<<<w.prop.multiProcessorCount*o.completion_blocks,128,0,w.stream>>>(cv,w.state,w.out,w.capacity);
   kit::end_epoch<<<1,1,0,w.stream>>>(w.state,w.capacity);ck(cudaPeekAtLastError());}
 };
 cudaGraph_t graph{};cudaGraphExec_t exec{};auto gs_start=Clock::now();
 if(o.graph){ck(cudaStreamBeginCapture(w.stream,cudaStreamCaptureModeThreadLocal));chain();ck(cudaStreamEndCapture(w.stream,&graph));ck(cudaGraphInstantiate(&exec,graph,0));}
 s.setup+=ms(gs_start);cudaEvent_t a{},b{};if constexpr(KIT_REUSE_EVENTS){a=w.timer_begin;b=w.timer_end;}else{ck(cudaEventCreate(&a));ck(cudaEventCreate(&b));}ck(cudaEventRecord(a,w.stream));
 if(o.graph)ck(cudaGraphLaunch(exec,w.stream));else chain();
 ck(cudaMemcpyAsync(w.hs,w.state,sizeof(*w.hs),cudaMemcpyDeviceToHost,w.stream));ck(cudaEventRecord(b,w.stream));ck(cudaEventSynchronize(b));float span;ck(cudaEventElapsedTime(&span,a,b));s.queue+=span;
 std::ostringstream diagnostic;diagnostic<<"{\"coarse\":"<<coarse<<",\"fine\":"<<fine<<",\"queue_ms\":"<<std::setprecision(12)<<span<<",\"last_count\":"<<w.hs->count<<",\"t\":"<<w.hs->t<<",\"q\":"<<w.hs->q<<",\"qc\":"<<w.hs->qc<<",\"qf\":"<<w.hs->qf<<",\"B\":"<<cover.B<<",\"D\":"<<cover.D<<",\"epochs\":"<<w.hs->epochs_done<<"}";if(!s.group_diagnostics.empty())s.group_diagnostics+=",";s.group_diagnostics+=diagnostic.str();
 if(w.hs->failed)std::cerr<<"GPU_REFUSAL p="<<p<<" coarse="<<coarse<<" fine="<<fine<<" failed="<<w.hs->failed<<" count="<<w.hs->count<<" capacity="<<w.capacity<<" t="<<w.hs->t<<" q="<<w.hs->q<<"\n";
 require(!w.hs->failed,"GPU seed/filter/overflow error: result refused");require(w.hs->epochs_done==epochs.size(),"incomplete GPU-prepared coverage");
 m[coarse]=w.hs->best_c>>32;m[fine]=w.hs->best_f>>32;pre[coarse]=w.hs->prefix_c;pre[fine]=w.hs->prefix_f;wit[coarse]=0xffffffffu-(u32)w.hs->best_c;wit[fine]=0xffffffffu-(u32)w.hs->best_f;
 if(exec)ck(cudaGraphExecDestroy(exec));if(graph)ck(cudaGraphDestroy(graph));if constexpr(!KIT_REUSE_EVENTS){ck(cudaEventDestroy(a));ck(cudaEventDestroy(b));}++s.groups;
}

static void output(u32 p,u32 mask,const Options&o,u32 iteration,const u32*m,const u32*pre,const u32*wi,double total,const Metrics&s){
 std::cout<<std::setprecision(12)<<"{\"type\":\"measurement\",\"status\":\"PASS\",\"p\":"<<p<<",\"mask\":"<<mask<<",\"scope\":\""<<o.scope<<"\",\"engine\":\""<<o.engine<<"\",\"iteration\":"<<iteration<<",\"memory\":\""<<o.memory<<"\",\"threads\":"<<o.threads<<",\"epoch\":"<<o.epoch<<",\"capacity\":"<<o.capacity<<",\"prep\":\""<<o.prep<<"\",\"pilot\":"<<o.pilot<<",\"target\":"<<(o.target?"true":"false")<<",\"graph\":"<<(o.graph?"true":"false")<<",\"total_ms\":"<<total<<",\"prepare_ms\":"<<s.prep<<",\"transfer_ms\":"<<s.transfer<<",\"gpu_prepare_ms\":"<<s.gpu_prepare<<",\"setup_ms\":"<<s.setup<<",\"queue_span_ms\":"<<s.queue<<",\"verify_ms\":"<<s.verify<<",\"groups\":"<<s.groups<<",\"fallbacks\":"<<s.fallbacks<<",\"profiles\":[";
 bool comma=false;for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2))){if(comma)std::cout<<',';comma=true;
  u32 K=std::max(m[r]+1,((p-1)/2)%r==0?2*pre[r]+2:pre[r]+2);
  std::cout<<"{\"r\":"<<r<<",\"max\":"<<m[r]<<",\"prefix\":"<<pre[r]<<",\"K\":"<<K<<",\"witness\":"<<wi[r]<<'}';}
 std::cout<<"],\"resources\":["<<s.resources<<"],\"chunk_count\":"<<(1u<<KIT_CHUNK_SHIFT)<<",\"completion_blocks_per_sm\":"<<o.completion_blocks<<",\"shared_min_resident\":"<<o.shared_min_resident<<",\"group_diagnostics\":["<<s.group_diagnostics<<"]}\n"<<std::flush;
}
int frozen_benchmark_main(int argc,char**argv){try{Options o=options(argc,argv);std::unique_ptr<Workspace>w;
 auto init=Clock::now();
 if(o.engine!="cpu"){ck(cudaFree(nullptr));cudaDeviceProp prop;ck(cudaGetDeviceProperties(&prop,0));require(prop.major==12&&prop.minor==0,"SM120 target required");
  std::cout<<"{\"type\":\"device\",\"sm\":120,\"multiprocessors\":"<<prop.multiProcessorCount<<",\"shared_optin\":"<<prop.sharedMemPerBlockOptin<<"}\n";
  if(o.engine=="candidate")w=std::make_unique<Workspace>(o.capacity);else baseline_init(*std::max_element(o.primes.begin(),o.primes.end()));}
 std::cout<<"{\"type\":\"initialization\",\"ms\":"<<ms(init)<<"}\n";
 if(o.profile)ck(cudaProfilerStart());
 for(unsigned it=0;it<o.repeats;++it)for(u32 p:o.primes){auto start=Clock::now();u32 mask=kit::requested_mask(p,o.scope),m[10]{},pre[10]{},wi[10]{};Metrics s;
  if(o.engine=="baseline"){auto a=Clock::now();require(!baseline_profile(p,mask,m,pre,o.scope=="all"),"baseline adapter mask");ck(cudaDeviceSynchronize());s.queue=ms(a);}
  else if(o.engine=="cpu"){kit::dense_host(p,mask,m,pre,wi);}
  else if(o.prep=="gpu"){auto cover=prepare_gpu_table(*w,p,mask,o,s);for(auto [c,f]:kit::groups(p,o.scope)){
   bool split=c<f&&((c==2&&(o.pair_mode&1))||(c==3&&(o.pair_mode&2)));
   if(split){gpu_prepared_group(*w,p,c,c,group_cover(p,c,c,cover,o),o,m,pre,wi,s);gpu_prepared_group(*w,p,f,f,group_cover(p,f,f,cover,o),o,m,pre,wi,s);}
   else gpu_prepared_group(*w,p,c,f,group_cover(p,c,f,cover,o),o,m,pre,wi,s);
  }}
  else for(auto [c,f]:kit::groups(p,o.scope)){
   // A genuine inability to seed small high-colour cases has a visible exact fallback.
   std::unique_ptr<kit::Prep> prep;
   try{prep=std::make_unique<kit::Prep>(p,c,f,o.pilot);}catch(const std::exception&e){
    if(std::string(e.what())!="fallback required"||p>2000000)throw;
    kit::dense_host(p,(1u<<(c-2))|(1u<<(f-2)),m,pre,wi);++s.fallbacks;continue;}
   s.prep+=prep->seconds*1000;run_group(*w,*prep,o,m,pre,wi,s);
  }
  auto verification=Clock::now();
  if(o.engine=="candidate")for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2))){AnchorStats z;z.maximum=m[r];z.prefix=pre[r];z.witness=wi[r];certify_witness(p,r,z,true);
   u32 precheck=1;while(precheck+1<p&&power(precheck+1,(p-1)/r,p)==1)++precheck;require(precheck==pre[r],"prefix direct verification");}
  if(p<=2000000){u32 dm[10]{},dp[10]{},dw[10]{};kit::dense_host(p,mask,dm,dp,dw);for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2))){require(m[r]==dm[r]&&pre[r]==dp[r],"GPU/dense CPU mismatch");if(o.engine=="candidate")require(wi[r]==dw[r],"minimal witness/dense");}}
  s.verify=ms(verification);output(p,mask,o,it,m,pre,wi,ms(start),s);
 }
 if(o.profile)ck(cudaProfilerStop());return 0;
 }catch(const std::exception&e){std::cerr<<"VDW_TEST_FAILURE: "<<e.what()<<'\n';return 2;}}
