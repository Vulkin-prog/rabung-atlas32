#define main atlas_unused_main
#include "stream.cpp"
#undef main
int main(int argc,char**argv){try{
 require(argc==5,"usage scan_dir frames windows output");
 auto began=std::chrono::steady_clock::now();std::ifstream wi(argv[3]);U l,h;std::vector<std::pair<U,U>> windows;while(wi>>l>>h)windows.emplace_back(l,h);
 std::vector<std::array<std::array<U,64>,8>> hh(windows.size());
 std::vector<U> np(windows.size());std::ifstream mf(argv[2]),sf;std::string name,an,current;U id,lo,hi,count,off,size,ao,as,blocks=0,n=0;
 size_t window=0;
 while(mf>>id>>lo>>hi>>count>>name>>off>>size>>an>>ao>>as){Bytes f=read(sf,argv[1],name,current,off,size);require(u64(f.data()+24)==lo && u64(f.data()+32)==hi && u32(f.data()+40)==count && u32(f.data()+44)==1,"header");Bytes raw=unpack(f,80,u64(f.data()+8),u64(f.data()+16),48);U p=lo-1;size_t at=0;
  for(U j=0;j<count;j++){p+=var(raw,at);unsigned mask=raw.at(at++);while(window<windows.size() && p>=windows[window].second)window++;bool take=window<windows.size() && p>=windows[window].first && p<windows[window].second;
   if(take){np[window]++;n++;}for(int i=0;i<8;i++)if(mask&(1<<i)){U m=var(raw,at);var(raw,at);var(raw,at);if(take)hh[window][i][std::min(U(63),m)]++;}
  }require(at==raw.size(),"trailing bytes");blocks++;
 }
 std::ofstream o(argv[4]);o<<"{\"primes\":"<<n<<",\"blocks\":"<<blocks<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<",\"windows\":[";
 for(size_t j=0;j<windows.size();j++){if(j)o<<',';o<<"{\"lo\":"<<windows[j].first<<",\"hi\":"<<windows[j].second<<",\"primes\":"<<np[j]<<",\"hist\":{";
  for(int i=0;i<8;i++){if(i)o<<',';o<<'"'<<i+2<<"\":[";for(int m=0;m<64;m++){if(m)o<<',';o<<hh[j][i][m];}o<<']';}o<<"}}";
 }o<<"]}\n";std::cerr<<"overlap complete primes="<<n<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
