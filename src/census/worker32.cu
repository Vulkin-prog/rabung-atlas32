// Frozen engine, new persistent input loop. Kernels and scientific functions are unchanged.
#include "engine32/src/frozen_entry.cu"
static void campaign_output(u64 p,u32 mask,const Options&o,u32 iteration,const u32*m,const u32*pre,const u32*wi,double total,const Metrics&s){
 std::cout<<std::setprecision(12)<<"{\"type\":\"measurement\",\"status\":\"PASS\",\"p\":"<<p<<",\"mask\":"<<mask<<",\"total_ms\":"<<total<<",\"profiles\":[";
 bool comma=false;for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2))){if(comma)std::cout<<',';comma=true;
  u64 K=std::max<u64>((u64)m[r]+1,((p-1)/2)%r==0?2ull*pre[r]+2:(u64)pre[r]+2);
  std::cout<<"{\"r\":"<<r<<",\"max\":"<<m[r]<<",\"prefix\":"<<pre[r]<<",\"K\":"<<K<<",\"witness\":"<<wi[r]<<'}';
 }std::cout<<"]}\n"<<std::flush;
}
int main(int argc,char**argv){try{
 Options o=options(argc,argv);require(o.engine=="candidate"&&o.prep=="gpu"&&o.repeats==1&&!o.profile,"production GPU options");
 ck(cudaFree(nullptr));Workspace w(o.capacity);
 std::cout<<"{\"type\":\"ready\",\"domain_bits\":32}\n"<<std::flush;
 std::string line;u64 expected=0;
 while(std::getline(std::cin,line)){
  if(line=="QUIT"){std::cout<<"{\"type\":\"bye\"}\n"<<std::flush;return 0;}
  std::istringstream in(line);std::string token;u64 id;require(bool(in>>id)&&id==expected,"batch sequence");std::vector<u32> primes;
  while(in>>token){require(token.find_first_not_of("0123456789")==std::string::npos,"prime token");u64 value=std::stoull(token);require(value>=1000&&value<(1ull<<32),"worker prime range");u32 p=(u32)value;require(prime(p),"worker primality");kit::requested_mask(p,o.scope);primes.push_back(p);}
  require(!primes.empty()&&primes.size()<=128,"worker batch size");auto batch_start=Clock::now();
  for(u32 p:primes){auto start=Clock::now();u32 mask=kit::requested_mask(p,o.scope),m[10]{},pre[10]{},wi[10]{};Metrics s;
  if(o.engine=="baseline"){auto a=Clock::now();require(!baseline_profile(p,mask,m,pre,o.scope=="all"),"baseline adapter mask");ck(cudaDeviceSynchronize());s.queue=ms(a);}
  else if(o.engine=="cpu"){kit::dense_host(p,mask,m,pre,wi);}
  else if(o.prep=="gpu"){auto cover=prepare_gpu_table(w,p,mask,o,s);for(auto [c,f]:kit::groups(p,o.scope)){
   bool split=c<f&&((c==2&&(o.pair_mode&1))||(c==3&&(o.pair_mode&2)));
   if(split){gpu_prepared_group(w,p,c,c,group_cover(p,c,c,cover,o),o,m,pre,wi,s);gpu_prepared_group(w,p,f,f,group_cover(p,f,f,cover,o),o,m,pre,wi,s);}
   else gpu_prepared_group(w,p,c,f,group_cover(p,c,f,cover,o),o,m,pre,wi,s);
  }}
  else for(auto [c,f]:kit::groups(p,o.scope)){
   // A genuine inability to seed small high-colour cases has a visible exact fallback.
   std::unique_ptr<kit::Prep> prep;
   try{prep=std::make_unique<kit::Prep>(p,c,f,o.pilot);}catch(const std::exception&e){
    if(std::string(e.what())!="fallback required"||p>2000000)throw;
    kit::dense_host(p,(1u<<(c-2))|(1u<<(f-2)),m,pre,wi);++s.fallbacks;continue;}
   s.prep+=prep->seconds*1000;run_group(w,*prep,o,m,pre,wi,s);
  }
  auto verification=Clock::now();
  if(o.engine=="candidate")for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2))){AnchorStats z;z.maximum=m[r];z.prefix=pre[r];z.witness=wi[r];certify_witness(p,r,z,true);
   u32 precheck=1;while(precheck+1<p&&power(precheck+1,(p-1)/r,p)==1)++precheck;require(precheck==pre[r],"prefix direct verification");}
 
  s.verify=ms(verification);campaign_output(p,mask,o,0,m,pre,wi,ms(start),s);
  }
  std::cout<<"{\"type\":\"batch_done\",\"id\":"<<id<<",\"count\":"<<primes.size()<<",\"wall_ms\":"<<ms(batch_start)<<"}\n"<<std::flush;++expected;
 }
 return 0;
 }catch(const std::exception&e){std::cerr<<"VDW_CAMPAIGN_FAILURE: "<<e.what()<<'\n';return 2;}}
