#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
// Independent full-integer segmented sieve, distinct from the Python odd sieve.
int main(int argc,char**argv){try{
 if(argc!=3)throw std::runtime_error("lo hi required");
 uint64_t lo=std::stoull(argv[1]),hi=std::stoull(argv[2]);
 if(lo<1||hi<lo||hi>(1ull<<33)||hi-lo>2000000)throw std::runtime_error("bounded sieve domain");
 uint64_t top=std::sqrt((long double)(hi?hi-1:0));while((top+1)*(top+1)<hi)++top;while(top*top>=hi&&top)--top;
 std::vector<unsigned char> base(top+1,1),mark(hi-lo,1);
 if(base.size())base[0]=0;if(base.size()>1)base[1]=0;
 for(uint64_t p=2;p<=top;++p)if(base[p]){
  if(p*p<=top)for(uint64_t j=p*p;j<=top;j+=p)base[j]=0;
  uint64_t first=std::max(p*p,((lo+p-1)/p)*p);
  for(uint64_t j=first;j<hi;j+=p)mark[j-lo]=0;
 }
 for(uint64_t p=lo;p<hi;++p)if(p>=2&&mark[p-lo])std::cout.write(reinterpret_cast<const char*>(&p),8);
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
