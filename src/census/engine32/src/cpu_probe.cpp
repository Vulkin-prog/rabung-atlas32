#include "preparation.hpp"
int main(int argc,char**argv){try{require(argc==4,"usage: cpu_probe prime scope pilot");u32 p=parse(argv[1]),pilot=parse(argv[3]);std::string scope=argv[2];
 auto a=std::chrono::steady_clock::now();double sec=0;unsigned groups=0;for(auto [c,f]:kit::groups(p,scope)){
  kit::Prep prep(p,c,f,pilot);sec+=prep.seconds;++groups;std::cout<<"{\"type\":\"prepare_group\",\"p\":"<<p<<",\"c\":"<<c<<",\"f\":"<<f<<",\"B\":"<<prep.B<<",\"D\":"<<prep.D<<",\"seconds\":"<<prep.seconds<<"}\n";
 }std::cout<<"{\"type\":\"prepare_total\",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-a).count()<<",\"groups\":"<<groups<<",\"gpu_executed\":false}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
