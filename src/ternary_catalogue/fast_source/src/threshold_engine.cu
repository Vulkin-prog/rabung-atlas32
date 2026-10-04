#include "gpu_pipeline.cuh"
#include "fast_capture.cuh"
#include "small_table.cuh"
#include <cuda_profiler_api.h>
#include <sstream>
#include <cstdio>
static void ck(cudaError_t e){if(e!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(e));}
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a){return std::chrono::duration<double,std::milli>(Clock::now()-a).count();}
struct Options {std::string engine="candidate",scope="all",memory="shared",prep="gpu";unsigned threshold_fixed_q=0,threshold_r=3,threshold_k=0,threshold_ab=0,threshold_singleton=0,grid_policy=3,cover_q=32,denominator_divisor=3,generator_blocks=4,completion_blocks=12,other_q=0,group36_q=0,single4_q=0,pair_mode=0,gpu_pilot=0,dynamic_t=0,t_mode=0,pilot_threads=0,shared_min_resident=0,threads=512,epoch=8192,pilot=64,repeats=1,capacity=2097152;bool graph=false,target=true,profile=false;std::vector<u64> primes;};
static Options options(int argc,char**argv){Options o;
 for(int i=1;i<argc;++i){std::string k=argv[i];if(k=="--profile"){o.profile=true;continue;}if(i+1==argc)throw std::runtime_error("missing option value");std::string v=argv[++i];
  if(k=="--threshold-fixed-q")o.threshold_fixed_q=parse(v.c_str());else if(k=="--threshold-r")o.threshold_r=parse(v.c_str());else if(k=="--threshold-k")o.threshold_k=parse(v.c_str());else if(k=="--threshold-ab")o.threshold_ab=parse(v.c_str());else if(k=="--threshold-singleton")o.threshold_singleton=parse(v.c_str());else if(k=="--other-q")o.other_q=parse(v.c_str());else if(k=="--group36-q")o.group36_q=parse(v.c_str());else if(k=="--single4-q")o.single4_q=parse(v.c_str());else if(k=="--pair-mode")o.pair_mode=parse(v.c_str());else if(k=="--t-mode")o.t_mode=parse(v.c_str());else if(k=="--pilot-threads")o.pilot_threads=parse(v.c_str());else if(k=="--gpu-pilot")o.gpu_pilot=parse(v.c_str());else if(k=="--dynamic-t")o.dynamic_t=parse(v.c_str());else if(k=="--grid-policy")o.grid_policy=parse(v.c_str());else if(k=="--cover-q")o.cover_q=parse(v.c_str());else if(k=="--denominator-divisor")o.denominator_divisor=parse(v.c_str());else if(k=="--generator-blocks")o.generator_blocks=parse(v.c_str());else if(k=="--completion-blocks")o.completion_blocks=parse(v.c_str());else if(k=="--shared-min-resident")o.shared_min_resident=parse(v.c_str());else if(k=="--engine")o.engine=v;else if(k=="--scope")o.scope=v;else if(k=="--memory")o.memory=v;else if(k=="--prep")o.prep=v;
  else if(k=="--threads")o.threads=parse(v.c_str());else if(k=="--epoch")o.epoch=parse(v.c_str());else if(k=="--pilot")o.pilot=parse(v.c_str());
  else if(k=="--repeat")o.repeats=parse(v.c_str());else if(k=="--capacity")o.capacity=parse(v.c_str());
  else if(k=="--graph"){require(v=="0"||v=="1","graph boolean");o.graph=v=="1";}
  else if(k=="--target"){require(v=="0"||v=="1","target boolean");o.target=v=="1";}
  else if(k=="--primes"){std::istringstream s(v);std::string q;while(std::getline(s,q,','))o.primes.push_back(parse64(q.c_str()));}
  else throw std::runtime_error("unknown option: "+k);
 }
 require((o.threshold_r==3||o.threshold_r==6)&&o.threshold_k<=256&&o.threshold_ab<=1&&o.threshold_singleton<=1,"threshold options");require(!o.threshold_k||(!o.graph&&o.scope=="36"),"threshold requires scope36 and graph0");
 require(o.threshold_fixed_q<=1&&(!o.threshold_fixed_q||(o.threshold_r==3&&o.threshold_singleton&&o.threshold_k>=4&&o.threshold_k-1<=(o.group36_q?o.group36_q:o.cover_q))),"fixed decision filter requires singleton3, k>=4 and coverQ>=k-1");
 require(!o.primes.empty()&&o.primes.size()<=128,"1..128 primes required");
 require(o.engine=="candidate","wide33 supports candidate GPU engine only");
 require(o.prep=="gpu","wide33 supports GPU preparation only");require(o.memory=="shared"||o.memory=="l2","memory");require(o.threads>=128&&o.threads<=1024&&o.threads%32==0,"threads");
 require(o.grid_policy<=3,"grid policy");
 require(o.gpu_pilot<=8192&&o.dynamic_t<=1,"bounded GPU pilot / dynamic threshold");
 require(o.pair_mode<=3,"pair splitting mode");
 require(!o.single4_q||(o.single4_q>=3&&o.single4_q<=o.cover_q&&(o.pair_mode&1)&&o.prep=="gpu"),"single4 cover requires GPU split24 and 3<=Q4<=Q");
 require(!o.group36_q||(o.group36_q>=3&&o.group36_q<=o.cover_q&&!(o.pair_mode&2)&&o.prep=="gpu"),"36 cover requires coupled GPU group and 3<=Q36<=Q");
 require(!o.other_q||(o.other_q>=3&&o.other_q<=o.cover_q&&o.prep=="gpu"),"singleton cover requires GPU prep and 3<=Qr<=Q");
 require(o.t_mode<=13&&(o.pilot_threads==0||o.pilot_threads==128||o.pilot_threads==256||o.pilot_threads==512||o.pilot_threads==1024),"anchor mode / pilot threads");
 require(o.cover_q>=3&&o.cover_q<=256&&o.denominator_divisor>=1&&o.denominator_divisor<=16,"cover options");
 require(o.prep=="gpu"||(o.cover_q==16&&o.denominator_divisor==2),"cover tuning is GPU-only");
 require(o.generator_blocks>=1&&o.generator_blocks<=8&&o.completion_blocks>=1&&o.completion_blocks<=64&&o.shared_min_resident<=3,"bounded experiments");
 require(o.epoch>=32&&o.epoch<=65536&&o.pilot<=512&&o.repeats>=1&&o.repeats<=16&&o.capacity>=1&&o.capacity<=4194304,"bounded options");
 for(u64 p:o.primes){require(p>=1000&&p<(1ull<<34)&&wide33::prime(p),"prime domain [1000,2^33)");kit::requested_mask(p,o.scope);}return o;
}
// Disjoint denominator intervals retain complete coverage exactly once.
static std::vector<std::pair<u32,u32>> gpu_epochs(u32 D,u32 epoch,u32 pilot){
 std::vector<std::pair<u32,u32>> out;u32 lo=1;
 if(pilot&&pilot<D){out.push_back({1,pilot+1});lo=pilot+1;}
 for(;lo<=D;lo+=epoch)out.push_back({lo,std::min(D+1,lo+epoch)});
 return out;
}
struct Workspace {
 cudaStream_t stream{};kit::State *state{},*hs{};kit::Candidate*out{};unsigned char*colour{},*hcolour{};u64 *inverse{},*hinverse{};
 u32 *cp[4]{},*fp[4]{},*hcp[4]{},*hfp[4]{};size_t plan_words=150000;unsigned capacity;cudaDeviceProp prop{};unsigned short*base{};unsigned*base_error{};u64*roots{},*hroots{};u32*labels{},*hlabels{};
 Workspace(unsigned cap):capacity(cap){ck(cudaGetDeviceProperties(&prop,0));require(prop.major==12&&prop.minor==0,"this kit is qualified only for SM120");
  ck(cudaMalloc(&base,1048577*2ull));ck(cudaMalloc(&base_error,4));ck(cudaMalloc(&roots,2560*8));ck(cudaMalloc(&labels,2560*4));ck(cudaMallocHost(&hroots,2560*8));ck(cudaMallocHost(&hlabels,2560*4));
  ck(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));ck(cudaMalloc(&state,sizeof(*state)));ck(cudaMallocHost(&hs,sizeof(*hs)));
  ck(cudaMalloc(&out,(size_t)cap*sizeof(*out)));ck(cudaMalloc(&colour,1048577));ck(cudaMallocHost(&hcolour,1048577));
  ck(cudaMalloc(&inverse,1048577*8ull));ck(cudaMallocHost(&hinverse,1048577*8ull));
  for(int i=0;i<4;++i){ck(cudaMalloc(&cp[i],plan_words*4));ck(cudaMalloc(&fp[i],plan_words*4));ck(cudaMallocHost(&hcp[i],plan_words*4));ck(cudaMallocHost(&hfp[i],plan_words*4));}
 }
 ~Workspace(){cudaStreamSynchronize(stream);cudaFree(base);cudaFree(base_error);cudaFree(roots);cudaFree(labels);cudaFreeHost(hroots);cudaFreeHost(hlabels);cudaStreamSynchronize(stream);for(int i=0;i<4;++i){cudaFree(cp[i]);cudaFree(fp[i]);cudaFreeHost(hcp[i]);cudaFreeHost(hfp[i]);}cudaFree(state);cudaFree(out);cudaFree(colour);cudaFree(inverse);cudaFreeHost(hs);cudaFreeHost(hcolour);cudaFreeHost(hinverse);cudaStreamDestroy(stream);}
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
struct Metrics{bool early_rejected=false;std::string group_diagnostics;double prep=0,transfer=0,gpu_prepare=0,setup=0,queue=0,verify=0;unsigned groups=0,fallbacks=0;std::string resources;};
static vdw_clipped::Cover prime_cover(u64 p,u32 Q,const Options&o){
 u64 area=(u64)(Q-1)*p,B=(u64)std::sqrt((long double)area);
 while((B+1)*(B+1)<=area)++B;while(B*B>area)--B;++B;
 u32 D=std::max(1u,(u32)(2*B/(Q-1))/o.denominator_divisor);
 while((u64)(Q-1)*D*D>4ull*p)--D;
 auto cover=vdw_clipped::from_denominator(p,Q,D);
 return cover;
}
static vdw_clipped::Cover prepare_gpu_table(Workspace&w,u64 p,u32 mask,const Options&o,Metrics&s){
 auto a=Clock::now();auto cover=prime_cover(p,o.cover_q,o);
 require(cover.B<=1048576,"wide34 table capacity");u32 inverse_D=cover.D;
 if(o.single4_q&&(mask&(1u<<(4-2)))){auto c4=prime_cover(p,o.single4_q,o);require(c4.B<=cover.B,"single4 table beyond prepared base");inverse_D=std::max(inverse_D,c4.D);}
 if(o.group36_q&&(mask&((1u<<(3-2))|(1u<<(6-2))))){auto c36=prime_cover(p,o.group36_q,o);require(c36.B<=cover.B,"36 table beyond prepared base");inverse_D=std::max(inverse_D,c36.D);}
 if(o.other_q&&(mask&(232u|((o.pair_mode&2)?16u:0u)))){auto co=prime_cover(p,o.other_q,o);require(co.B<=cover.B,"singleton table beyond prepared base");inverse_D=std::max(inverse_D,co.D);}
 require(inverse_D<=1048576,"wide34 inverse capacity");u32 m=1;for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2)))m=std::lcm(m,r);
 u64 ni=wide33::ninverse(p),one=((unsigned __int128)1<<wide33::radix_bits)%p,R2=(unsigned __int128)one*one%p;
 #if KIT_CHARACTER_ROOT
 u64 z=character_root(p,m);
#else
 u64 z=wide33::power128(root(p),(p-1)/m,p);
#endif
 u64 cur=1;std::vector<std::pair<u64,u32>> v(m);
 for(u32 j=0;j<m;++j){v[j]={wide33::mont(cur,R2,p,ni),j};cur=(unsigned __int128)cur*z%p;}std::sort(v.begin(),v.end());
 for(u32 j=0;j<m;++j){w.hroots[j]=v[j].first;w.hlabels[j]=v[j].second;}
 s.prep+=ms(a);cudaEvent_t b,e;ck(cudaEventCreate(&b));ck(cudaEventCreate(&e));ck(cudaEventRecord(b,w.stream));
 ck(cudaMemcpyAsync(w.roots,w.hroots,m*8,cudaMemcpyHostToDevice,w.stream));ck(cudaMemcpyAsync(w.labels,w.hlabels,m*4,cudaMemcpyHostToDevice,w.stream));ck(cudaMemsetAsync(w.base_error,0,4,w.stream));
 kit::small_table<<<(cover.B+256)/256,256,0,w.stream>>>(p,cover.B,m,ni,R2,one,w.roots,w.labels,w.base,w.base_error);
 kit::inverse_table<<<(inverse_D+256)/256,256,0,w.stream>>>(p,inverse_D,ni,R2,one,w.inverse);ck(cudaPeekAtLastError());
 ck(cudaEventRecord(e,w.stream));ck(cudaEventSynchronize(e));float dt;ck(cudaEventElapsedTime(&dt,b,e));s.gpu_prepare+=dt;
 ck(cudaEventDestroy(b));ck(cudaEventDestroy(e));return cover;
}
static void gpu_prepared_group(Workspace&w,u64 p,u32 coarse,u32 fine,vdw_clipped::Cover cover,const Options&o,u32*m,u32*pre,u64*wit,Metrics&s){
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
 w.hs->stop_k=0; // Catalogue capture must never stop early.w.hs->stop_r=o.threshold_r;w.hs->stop_ab=o.threshold_ab;w.hs->stop_fusion=((p-1)/2)%o.threshold_r==0;w.hs->fixed_q=o.threshold_fixed_q?o.threshold_k-1:0;
 vdw_next::CompletionView cv{p,fine,cover.B,3,cover.D,w.colour,w.inverse};
#if VDW_COMPLETION_MONT
 cv.ni=wide33::ninverse(p);
#endif
 cv.sigma_value=((p-1)/2)%fine;
auto epochs=gpu_epochs(cover.D,o.epoch,o.gpu_pilot);
s.setup+=ms(host_start);
 auto threshold_reached=[&](bool after_pilot){
  if(o.threshold_k){kit::threshold_stop<<<1,1,0,w.stream>>>(w.state);ck(cudaPeekAtLastError());}
  // Only the fixed18 large-prime path observes the flag after the pilot.
  // No speculative result: the preceding completed same-stream epoch proves it.
  if(after_pilot&&o.threshold_fixed_q&&o.threshold_k==18&&p>=(1ull<<32)){
   ck(cudaMemcpyAsync(w.hs,w.state,sizeof(*w.hs),cudaMemcpyDeviceToHost,w.stream));
   ck(cudaStreamSynchronize(w.stream));
   return w.hs->stop_flag!=0;
  }
  return false;
 };
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
  u32 segments=(cover.B+1023)/1024;kit::seed_runs<<<KIT_WARP_SEED?(segments+3)/4:(segments+127)/128,128,0,w.stream>>>(p,cover.B,coarse,fine,w.colour,w.state);
  kit::finish_seed<<<1,1,0,w.stream>>>(cv,w.state,w.base_error,o.t_mode);
  if(threshold_reached(false))return;
  for(auto [lo,hi]:epochs){
   kit::epoch_snapshot<<<1,1,0,w.stream>>>(w.state,lo,hi,o.dynamic_t&&lo>1,o.t_mode);
   Options launch=o;if(lo==1&&o.gpu_pilot&&o.pilot_threads)launch.threads=o.pilot_threads;
   for(int i=0;i<4;++i)ck(kit::launch_group(gs[i],c[i],f[i],w.state,selected,o.target,w.out,w.capacity,grid_for_epoch(launch,w,cover.D,lo,hi),launch.threads,shared[i],w.stream));
   kit::capture_candidates<<<w.prop.multiProcessorCount*o.completion_blocks,128,0,w.stream>>>(cv,w.state,w.out,w.capacity);
   kit::end_epoch<<<1,1,0,w.stream>>>(w.state,w.capacity);ck(cudaPeekAtLastError());if(threshold_reached(lo==1))return;}
 };
 cudaGraph_t graph{};cudaGraphExec_t exec{};auto gs_start=Clock::now();
 if(o.graph){ck(cudaStreamBeginCapture(w.stream,cudaStreamCaptureModeThreadLocal));chain();ck(cudaStreamEndCapture(w.stream,&graph));ck(cudaGraphInstantiate(&exec,graph,0));}
 s.setup+=ms(gs_start);cudaEvent_t a,b;ck(cudaEventCreate(&a));ck(cudaEventCreate(&b));ck(cudaEventRecord(a,w.stream));
 if(o.graph)ck(cudaGraphLaunch(exec,w.stream));else chain();
 ck(cudaMemcpyAsync(w.hs,w.state,sizeof(*w.hs),cudaMemcpyDeviceToHost,w.stream));ck(cudaEventRecord(b,w.stream));ck(cudaEventSynchronize(b));float span;ck(cudaEventElapsedTime(&span,a,b));s.queue+=span;s.early_rejected=w.hs->stop_flag;
 std::ostringstream diagnostic;diagnostic<<"{\"coarse\":"<<coarse<<",\"fine\":"<<fine<<",\"queue_ms\":"<<std::setprecision(12)<<span<<",\"last_count\":"<<w.hs->count<<",\"t\":"<<w.hs->t<<",\"q\":"<<w.hs->q<<",\"qc\":"<<w.hs->qc<<",\"qf\":"<<w.hs->qf<<",\"B\":"<<cover.B<<",\"D\":"<<cover.D<<",\"epochs\":"<<w.hs->epochs_done<<"}";if(!s.group_diagnostics.empty())s.group_diagnostics+=",";s.group_diagnostics+=diagnostic.str();
 if(w.hs->failed)std::cerr<<"GPU_REFUSAL p="<<p<<" coarse="<<coarse<<" fine="<<fine<<" failed="<<w.hs->failed<<" count="<<w.hs->count<<" capacity="<<w.capacity<<" t="<<w.hs->t<<" q="<<w.hs->q<<"\n";
 require(!w.hs->failed,"GPU seed/filter/overflow error: result refused");require(s.early_rejected||w.hs->epochs_done==epochs.size(),"incomplete GPU-prepared coverage");
 m[coarse]=w.hs->best_c>>48;m[fine]=w.hs->best_f>>48;pre[coarse]=w.hs->prefix_c;pre[fine]=w.hs->prefix_f;wit[coarse]=0xffffffffffffull-(w.hs->best_c&0xffffffffffffull);wit[fine]=0xffffffffffffull-(w.hs->best_f&0xffffffffffffull);
 if(exec)ck(cudaGraphExecDestroy(exec));if(graph)ck(cudaGraphDestroy(graph));ck(cudaEventDestroy(a));ck(cudaEventDestroy(b));++s.groups;
}

static void output(u64 p,u32 mask,const Options&o,u32 iteration,const u32*m,const u32*pre,const u64*wi,double total,const Metrics&s){
 std::cout<<std::setprecision(12)<<"{\"type\":\"measurement\",\"status\":\"PASS\",\"p\":"<<p<<",\"mask\":"<<mask<<",\"scope\":\""<<o.scope<<"\",\"engine\":\""<<o.engine<<"\",\"iteration\":"<<iteration<<",\"memory\":\""<<o.memory<<"\",\"threads\":"<<o.threads<<",\"epoch\":"<<o.epoch<<",\"capacity\":"<<o.capacity<<",\"prep\":\""<<o.prep<<"\",\"pilot\":"<<o.pilot<<",\"target\":"<<(o.target?"true":"false")<<",\"graph\":"<<(o.graph?"true":"false")<<",\"total_ms\":"<<total<<",\"prepare_ms\":"<<s.prep<<",\"transfer_ms\":"<<s.transfer<<",\"gpu_prepare_ms\":"<<s.gpu_prepare<<",\"setup_ms\":"<<s.setup<<",\"queue_span_ms\":"<<s.queue<<",\"verify_ms\":"<<s.verify<<",\"groups\":"<<s.groups<<",\"fallbacks\":"<<s.fallbacks<<",\"profiles\":[";
 bool comma=false;for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2))){if(comma)std::cout<<',';comma=true;
  u64 K=std::max<u64>((u64)m[r]+1,((p-1)/2)%r==0?2ull*pre[r]+2:(u64)pre[r]+2);
  std::cout<<"{\"r\":"<<r<<",\"max\":"<<m[r]<<",\"prefix\":"<<pre[r]<<",\"K\":"<<K<<",\"witness\":"<<wi[r]<<'}';}
 std::cout<<"],\"resources\":["<<s.resources<<"],\"chunk_count\":"<<(1u<<KIT_CHUNK_SHIFT)<<",\"completion_blocks_per_sm\":"<<o.completion_blocks<<",\"shared_min_resident\":"<<o.shared_min_resident<<",\"group_diagnostics\":["<<s.group_diagnostics<<"]}\n"<<std::flush;
}
