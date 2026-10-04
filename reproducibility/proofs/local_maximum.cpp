#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using U=uint64_t;
int m,rank;std::vector<unsigned> column;std::vector<U> dual;
U syndrome(int k,bool iid=false){
 if(k<=1)return 0;
 unsigned ns=iid?1:(1u<<rank);int cap=k-1;
 std::vector<U>a(ns*2*cap),b(a.size());
 auto at=[&](unsigned s,int last,int run){return (s*2+last)*cap+run-1;};
 a[at(0,0,1)]=1;a[at(iid?0:column[0],1,1)]=1;
 for(int j=1;j<m;j++){
  std::fill(b.begin(),b.end(),0);unsigned mask=iid?0:column[j];
  for(unsigned s=0;s<ns;s++)for(int bit=0;bit<2;bit++)for(int r=1;r<=cap;r++){
   U v=a[at(s,bit,r)];if(!v)continue;
   b[at(s^((1-bit)?mask:0),1-bit,1)]+=v;
   if(r<cap)b[at(s^(bit?mask:0),bit,r+1)]+=v;
  }a.swap(b);
 }
 U total=0;for(int bit=0;bit<2;bit++)for(int r=1;r<=cap;r++)total+=a[at(0,bit,r)];return total;
}
U fourier(int k){
 if(k<=1)return 0;int cap=k-1;__int128 all=0;
 for(U u:dual){
  std::vector<int64_t>a(2*cap),b(2*cap);a[0]=1;a[cap]=(u&1)?-1:1;
  for(int j=1;j<m;j++){
   std::fill(b.begin(),b.end(),0);int sign=(u>>j&1)?-1:1;
   for(int bit=0;bit<2;bit++)for(int r=0;r<cap;r++){
    int64_t v=a[bit*cap+r];b[(1-bit)*cap]+=v*(bit?1:sign);
    if(r+1<cap)b[bit*cap+r+1]+=v*(bit?sign:1);
   }a.swap(b);
  }for(auto v:a)all+=v;
 }
 if(all<0||all%dual.size())throw std::runtime_error("Fourier divisibility");return U(all/dual.size());
}
int main(int argc,char**argv){try{
 if(argc<3)throw std::runtime_error("usage relation-input csv-output [one-k]");
 std::ifstream in(argv[1]);in>>m>>rank;std::vector<U>R(rank);column.assign(m,0);
 for(int i=0;i<rank;i++){in>>R[i];for(int j=0;j<m;j++)if(R[i]>>j&1)column[j]^=1u<<i;}
 if(!in)throw std::runtime_error("input");dual={0};for(U r:R){size_t n=dual.size();for(size_t i=0;i<n;i++)dual.push_back(dual[i]^r);}
 int low=argc==4?std::stoi(argv[3]):1,high=argc==4?low:m+1;
 std::ofstream out(argv[2]);out<<"threshold,below_code_syndrome,below_code_fourier,below_iid\n";
 auto t=std::chrono::steady_clock::now();
 for(int k=low;k<=high;k++){
  U a=syndrome(k),b=fourier(k),c=syndrome(k,true);if(a!=b)throw std::runtime_error("algorithm disagreement");
  out<<k<<','<<a<<','<<b<<','<<c<<'\n';
 }
 std::cout<<"PASS "<<high-low+1<<" thresholds two algorithms seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<'\n';
 }catch(std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
