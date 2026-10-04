#define main atlas_unused_main
#include "stream.cpp"
#undef main
int main(int argc,char**argv){try{
 require(argc==4,"usage scan_dir frames output");
 U counts[8][128][5]{};auto began=std::chrono::steady_clock::now();std::ifstream mf(argv[2]),sf;std::string name,an,current;U id,lo,hi,count,off,size,ao,as,nb=0;
 while(mf>>id>>lo>>hi>>count>>name>>off>>size>>an>>ao>>as){Bytes f=read(sf,argv[1],name,current,off,size);Bytes raw=unpack(f,80,u64(f.data()+8),u64(f.data()+16),48);U p=lo-1;size_t at=0;
  for(U row=0;row<count;row++){p+=var(raw,at);unsigned mask=raw.at(at++);
   for(int i=0;i<8;i++)if(mask&(1<<i)){U m=var(raw,at),l=var(raw,at),w=var(raw,at);if(m<U(cuts[i]))continue;require(m<128,"M capacity");U c=0;bool sigma=((p-1)/2)%(i+2)==0;
    if(sigma){if(l>1)c=l;else for(int q:small_primes){require(U(q)<=m+1,"center cap");if(power(q,(p-1)/(i+2),p)!=1){c=q-1;break;}}require(c && c<=m,"center bounds");}
    auto&v=counts[i][m];v[0]++;v[1]+=c==m;v[2]+=l==m;v[3]+=c==m || l==m;v[4]+=c==m && w==(p-c+1)/2;
   }
  }require(at==raw.size(),"trailing bytes");nb++;
 }
 std::ofstream o(argv[3]);o<<"{\"blocks\":"<<nb<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<",\"by_M\":{";
 for(int i=0;i<8;i++){if(i)o<<',';o<<'"'<<i+2<<"\":{";bool comma=false;for(int m=0;m<128;m++)if(counts[i][m][0]){if(comma)o<<',';comma=true;o<<'"'<<m<<"\":[";for(int j=0;j<5;j++){if(j)o<<',';o<<counts[i][m][j];}o<<']';}o<<'}';}o<<"}}\n";
 std::cerr<<"tails complete blocks="<<nb<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
