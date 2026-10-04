#define VDW_PERIODIC_NO_MAIN
#include "../src/periodic_reference.cpp"
int main(int argc,char**argv){try{u64 profiles=0,fallback=0,words=0,stages=0;
 bool medium=argc==2&&std::string(argv[1])=="--million";
 std::vector<u32> primes;
 if(medium){for(u32 p=1000000;primes.size()<8;++p)if(prime(p))primes.push_back(p);}
 else for(u32 p=11;p<450;p+=2)if(prime(p))primes.push_back(p);
 for(u32 p:primes)for(u32 r=2;r<=9;++r)if((p-1)%r==0)for(u32 div:{1u,2u,4u}){
  try{auto E=periodic_experiment(p,r,medium?64:16,true,div);require(E.exp.dense_checked,"expected dense check");++profiles;words+=E.exp.main.words;stages+=E.exp.main.stages;}
  catch(const std::invalid_argument&e){std::string msg=e.what();if(msg.find("fallback")==std::string::npos)throw;++fallback;}
  catch(const std::runtime_error&e){std::string msg=e.what();if(msg.find("fallback")==std::string::npos&&msg!="cover")throw;++fallback;}
 }
 std::cout<<"{\"status\":\"PASS\",\"domain\":\""<<(medium?"million":"small")<<"\",\"dense_profiles_compared\":"<<profiles<<",\"explicit_fallbacks\":"<<fallback<<",\"filter_words\":"<<words<<",\"filter_stages\":"<<stages<<",\"GPU_executed\":false}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
