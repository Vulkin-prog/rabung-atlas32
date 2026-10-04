#include "src/threshold_engine.cu"
int main(int argc,char**argv){try{
 u32 minimum=15;require(argc==3&&std::string(argv[1])=="--minimum","use --minimum value");minimum=parse(argv[2]);
 require(minimum>=4&&minimum<=24,"explicit retained threshold4..24");
 Options o;o.scope="36";o.threshold_singleton=1;o.threshold_fixed_q=1;o.threshold_k=minimum;o.threshold_r=3;
 o.cover_q=o.group36_q=o.other_q=minimum-1;o.denominator_divisor=3;o.memory="shared";o.threads=1024;o.epoch=16384;
 o.gpu_pilot=1536;o.pilot_threads=512;o.grid_policy=3;o.dynamic_t=0;o.t_mode=13;o.pair_mode=1;
 ck(cudaFree(nullptr));Workspace w(o.capacity);
 const u32 capacity=1048576;kit::CapturedRun*output{};unsigned long long*count{};
 ck(cudaMalloc(&output,(size_t)capacity*sizeof(kit::CapturedRun)));ck(cudaMalloc(&count,8));
 ck(cudaMemcpyToSymbol(kit::capture_output,&output,sizeof(output)));ck(cudaMemcpyToSymbol(kit::capture_count,&count,sizeof(count)));
 ck(cudaMemcpyToSymbol(kit::capture_capacity,&capacity,4));ck(cudaMemcpyToSymbol(kit::capture_minimum,&minimum,4));
 std::vector<kit::CapturedRun>runs(capacity);std::cout<<"{\"type\":\"ready\",\"minimum\":"<<minimum<<",\"domain_bits\":34}\n"<<std::flush;
 std::string line;u64 sequence=0;
 while(std::getline(std::cin,line)){
  if(line=="QUIT")break;std::istringstream in(line);u64 id;require(bool(in>>id)&&id==sequence,"batch sequence");
  std::vector<u64>primes;std::string token;while(in>>token){u64 p=parse64(token.c_str());require(p>=1000&&p<(1ull<<34)&&p%6==1&&wide33::prime(p),"prime domain");primes.push_back(p);}
  require(!primes.empty()&&primes.size()<=128,"bounded batch");auto batch_start=Clock::now();
  for(u64 p:primes){auto began=Clock::now();u32 maxima[10]{},prefix[10]{};u64 witness[10]{};Metrics metric;
   ck(cudaMemsetAsync(count,0,8,w.stream));auto cover=prepare_gpu_table(w,p,1u<<(3-2),o,metric);
   gpu_prepared_group(w,p,3,3,cover,o,maxima,prefix,witness,metric);
   require(!metric.early_rejected&&!w.hs->failed,"partial capture refused");unsigned long long n;
   ck(cudaMemcpy(&n,count,8,cudaMemcpyDeviceToHost));require(n<=capacity,"capture capacity");
   ck(cudaMemcpy(runs.data(),output,n*sizeof(kit::CapturedRun),cudaMemcpyDeviceToHost));
   std::sort(runs.begin(),runs.begin()+n,[](auto a,auto b){return a.start<b.start||(a.start==b.start&&a.length<b.length);});
   auto end=std::unique(runs.begin(),runs.begin()+n,[](auto a,auto b){require(a.start!=b.start||(a.length==b.length&&a.colour==b.colour),"duplicate capture disagreement");return a.start==b.start;});
   size_t unique=end-runs.begin();u64 ni=wide33::ninverse(p),natural_root=0;
   for(int j=0;j<3;++j)if(w.hlabels[j]==1)natural_root=wide33::mont(w.hroots[j],1,p,ni);
   require(natural_root&&natural_root!=1,"natural cubic root");
   for(size_t j=0;j<unique;++j){auto r=runs[j];require(r.start>=1&&r.length>=minimum&&r.start+r.length<p,"capture interval");}
   std::cout<<std::setprecision(12)<<"{\"type\":\"catalogue\",\"status\":\"PASS_COMPLETE_LONG_RUN_CATALOGUE\",\"p\":"<<p<<",\"r\":3,\"minimum\":"<<minimum<<",\"M_bound\":"<<maxima[3]<<",\"L\":"<<prefix[3]<<",\"complete_long_runs\":true,\"character_root\":"<<natural_root<<",\"raw_anchor_captures\":"<<n<<",\"unique_runs\":"<<unique<<",\"total_ms\":"<<ms(began)<<",\"gpu_prepare_ms\":"<<metric.gpu_prepare<<",\"queue_ms\":"<<metric.queue<<",\"B\":"<<cover.B<<",\"D\":"<<cover.D<<",\"runs\":[";
   for(size_t j=0;j<unique;++j){if(j)std::cout<<',';auto r=runs[j];std::cout<<'['<<r.start<<','<<r.length<<','<<r.colour<<']';}std::cout<<"]}\n"<<std::flush;
  }
  std::cout<<"{\"type\":\"batch_done\",\"id\":"<<id<<",\"count\":"<<primes.size()<<",\"wall_ms\":"<<ms(batch_start)<<"}\n"<<std::flush;++sequence;
 }
 cudaFree(output);cudaFree(count);return 0;
}catch(const std::exception&e){std::cerr<<"FAST_CATALOGUE_REFUSAL: "<<e.what()<<'\n';return 2;}}
