#include "src/threshold_engine.cu"
static void verify_bound33(u64 p,u32 r,u32 n,u32 prefix,u64 witness){
 require(n&&witness&&witness<=(p-1)/2&&(u64)witness+n<=p,"partial witness interval");u64 e=(p-1)/r,y=wide33::power128(witness,e,p);
 for(u32 j=1;j<n;++j)require(wide33::power128((u64)witness+j,e,p)==y,"partial witness colour");
 u32 l=1;while((u64)l+1<p&&wide33::power128((u64)l+1,e,p)==1)++l;require(l==prefix,"partial prefix verification");
}
int main(int argc,char**argv){try{
 Options o=options(argc,argv);require(o.repeats==1&&!o.profile,"persistent threshold options");
 ck(cudaFree(nullptr));Workspace w(o.capacity);std::cout<<"{\"type\":\"ready\",\"domain_bits\":34}\n"<<std::flush;
 std::string line;u64 expected=0;
 while(std::getline(std::cin,line)){
  if(line=="QUIT"){std::cout<<"{\"type\":\"bye\"}\n"<<std::flush;return 0;}
  std::istringstream in(line);u64 id;std::string token;require(bool(in>>id)&&id==expected,"batch sequence");std::vector<u64>primes;
  while(in>>token){u64 p=parse64(token.c_str());require(p>=1000&&p<(1ull<<34)&&wide33::prime(p)&&(p-1)%6==0,"threshold prime domain");primes.push_back(p);}
  require(!primes.empty()&&primes.size()<=128,"threshold batch size");auto bs=Clock::now();
  for(u64 p:primes){auto start=Clock::now();u32 m[10]{},pre[10]{};u64 wi[10]{};Metrics s;
   u32 mask=o.threshold_singleton?1u<<(o.threshold_r-2):kit::requested_mask(p,"36");
   u32 coarse=o.threshold_singleton?o.threshold_r:3,fine=o.threshold_singleton?o.threshold_r:6;
   // Very small domains can have M6<3, below the GPU anchor precondition.
   // Use an explicit complete CPU path there, rather than lowering an invariant.
   if(p<=2000000){u32 smallw[10]{};kit::dense_host(p,mask,m,pre,smallw);for(u32 j=0;j<10;++j)wi[j]=smallw[j];}
   else{auto cover=prepare_gpu_table(w,p,mask,o,s);gpu_prepared_group(w,p,coarse,fine,o.group36_q?prime_cover(p,o.group36_q,o):cover,o,m,pre,wi,s);}
   if(o.threshold_fixed_q&&p>2000000&&!s.early_rejected){
    // The fixed pass proves absence of k, but does not determine maxima below
    // k. Reuse the prepared tables for a separate original full-profile pass.
    Options full=o;full.threshold_fixed_q=0;full.threshold_k=0;full.threshold_ab=0;
    auto cover=prime_cover(p,o.group36_q?o.group36_q:o.cover_q,o);
    gpu_prepared_group(w,p,coarse,fine,cover,full,m,pre,wi,s);
    require(m[o.threshold_r]<o.threshold_k,"fixed threshold missed a run");
   }
   if(s.early_rejected)verify_bound33(p,o.threshold_r,m[o.threshold_r],pre[o.threshold_r],wi[o.threshold_r]);
   else{certify_witness33(p,coarse,m[coarse],pre[coarse],wi[coarse]);if(fine!=coarse)certify_witness33(p,fine,m[fine],pre[fine],wi[fine]);}
   u32 r=o.threshold_r;u64 border=((p-1)/2)%r==0?2ull*pre[r]+2:(u64)pre[r]+2;
   u64 boundK=std::max<u64>((u64)m[r]+1,border);bool accepted=o.threshold_k&&(o.threshold_ab?boundK<=o.threshold_k:m[r]<o.threshold_k);
   require(!accepted||!s.early_rejected,"partial result cannot be accepted");
   std::cout<<std::setprecision(12)<<"{\"type\":\"threshold\",\"status\":\""<<(s.early_rejected?"REJECT_BOUND":"COMPLETE")<<"\",\"engine\":\""<<(p<=2000000?"CPU_COMPLETE":"GPU")<<"\",\"p\":"<<p<<",\"r\":"<<r<<",\"k\":"<<o.threshold_k<<",\"ab\":"<<(o.threshold_ab?"true":"false")<<",\"accepted\":"<<(accepted?"true":"false")<<",\"M_bound\":"<<m[r]<<",\"L\":"<<pre[r]<<",\"K_bound\":"<<boundK<<",\"witness\":"<<wi[r]<<",\"complete\":"<<(!s.early_rejected?"true":"false")<<",\"total_ms\":"<<ms(start)<<",\"gpu_prepare_ms\":"<<s.gpu_prepare<<",\"queue_ms\":"<<s.queue<<",\"groups\":["<<s.group_diagnostics<<"]}\n"<<std::flush;
  }
  std::cout<<"{\"type\":\"batch_done\",\"id\":"<<id<<",\"count\":"<<primes.size()<<",\"wall_ms\":"<<ms(bs)<<"}\n"<<std::flush;++expected;
 }
 return 0;
}catch(const std::exception&e){std::cerr<<"THRESHOLD_REFUSAL: "<<e.what()<<'\n';return 2;}}
