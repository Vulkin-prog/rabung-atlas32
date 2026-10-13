#include <thread>
#include <atomic>
#include <mutex>
// Frozen engine, new persistent input loop. Kernels and scientific functions are unchanged.
#include "engine33/src/frozen_entry.cu"
static std::string campaign_output(u64 p,u32 mask,const Options&o,u32 iteration,const u32*m,const u32*pre,const u32*wi,double total,const Metrics&s){
 std::ostringstream result;result<<std::setprecision(12)<<"{\"type\":\"measurement\",\"status\":\"PASS\",\"p\":"<<p<<",\"mask\":"<<mask<<",\"total_ms\":"<<total<<",\"prepare_ms\":"<<s.prep<<",\"gpu_prepare_ms\":"<<s.gpu_prepare<<",\"setup_ms\":"<<s.setup<<",\"queue_ms\":"<<s.queue<<",\"verify_ms\":"<<s.verify<<",\"groups\":"<<s.groups<<",\"profiles\":[";
 bool comma=false;for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2))){if(comma)result<<',';comma=true;
  u64 K=std::max<u64>((u64)m[r]+1,((p-1)/2)%r==0?2ull*pre[r]+2:(u64)pre[r]+2);
  result<<"{\"r\":"<<r<<",\"max\":"<<m[r]<<",\"prefix\":"<<pre[r]<<",\"K\":"<<K<<",\"witness\":"<<wi[r]<<'}';
 }result<<"]}\n";return result.str();
}
static std::string calculate_prime(u64 p,const Options&o,Workspace&w){auto start=Clock::now();u32 mask=kit::requested_mask(p,o.scope),m[10]{},pre[10]{},wi[10]{};Metrics s;
 auto cover=prepare_gpu_table(w,p,mask,o,s);for(auto [c,f]:kit::groups(p,o.scope)){
   bool split=c<f&&((c==2&&(o.pair_mode&1))||(c==3&&(o.pair_mode&2)));
   if(split){gpu_prepared_group(w,p,c,c,cover,o,m,pre,wi,s);gpu_prepared_group(w,p,f,f,(f==4&&o.single4_q)?prime_cover(p,o.single4_q,o):(f>=5&&o.other_q)?prime_cover(p,o.other_q,o):cover,o,m,pre,wi,s);}
   else gpu_prepared_group(w,p,c,f,(c==4&&f==8&&o.single4_q)?prime_cover(p,o.single4_q,o):(c==3&&f==6&&o.group36_q)?prime_cover(p,o.group36_q,o):(c==f&&f>=5&&o.other_q)?prime_cover(p,o.other_q,o):cover,o,m,pre,wi,s);
  }
 auto verification=Clock::now();for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2)))certify_witness33(p,r,m[r],pre[r],wi[r]);

 s.verify=ms(verification);return campaign_output(p,mask,o,0,m,pre,wi,ms(start),s);
  }
int main(int argc,char**argv){try{
 Options o=options(argc,argv);require(o.engine=="candidate"&&o.prep=="gpu"&&o.repeats==1&&!o.profile,"production GPU options");
 ck(cudaFree(nullptr));std::vector<std::unique_ptr<Workspace>> pool;for(int k=0;k<4;++k)pool.push_back(std::make_unique<Workspace>(o.capacity));
 std::cout<<"{\"type\":\"ready\",\"domain_bits\":33}\n"<<std::flush;
 std::string line;u64 expected=0;
 while(std::getline(std::cin,line)){
  if(line=="QUIT"){std::cout<<"{\"type\":\"bye\"}\n"<<std::flush;return 0;}
  std::istringstream in(line);std::string token;u64 id;require(bool(in>>id)&&id==expected,"batch sequence");std::vector<u64> primes;
  while(in>>token){require(token.find_first_not_of("0123456789")==std::string::npos,"prime token");u64 value=std::stoull(token);require(value>=1000&&value<(1ull<<33),"worker prime range");u64 p=(u64)value;require(wide33::prime(p),"worker primality");kit::requested_mask(p,o.scope);primes.push_back(p);}
  require(!primes.empty()&&primes.size()<=128,"worker batch size");auto batch_start=Clock::now();

  std::vector<std::string> results(primes.size());std::atomic<size_t> next{0};std::atomic<bool> failed{false};std::mutex errors;std::string error;
  auto job=[&](size_t index){try{ck(cudaSetDevice(0));while(!failed){size_t k=next.fetch_add(1);if(k>=primes.size())break;results[k]=calculate_prime(primes[k],o,*pool[index]);}}
   catch(const std::exception&e){std::lock_guard<std::mutex> hold(errors);if(error.empty())error=e.what();failed=true;}};
  if constexpr(4==1)job(0);else{std::vector<std::thread> threads;for(size_t k=0;k<4;++k)threads.emplace_back(job,k);for(auto &t:threads)t.join();}
  if(failed)throw std::runtime_error(error);for(const auto &line:results)std::cout<<line;std::cout<<std::flush;
  std::cout<<"{\"type\":\"batch_done\",\"id\":"<<id<<",\"count\":"<<primes.size()<<",\"wall_ms\":"<<ms(batch_start)<<"}\n"<<std::flush;++expected;
 }
 return 0;
 }catch(const std::exception&e){std::cerr<<"VDW_CAMPAIGN_FAILURE: "<<e.what()<<'\n';return 2;}}
