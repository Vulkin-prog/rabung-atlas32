#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <openssl/sha.h>
#include <zlib.h>
using U=uint64_t; using I=int64_t; using Bytes=std::vector<unsigned char>;
U u64(const unsigned char* p){U v=0;for(int i=0;i<8;i++)v|=U(p[i])<<(8*i);return v;}
unsigned u32(const unsigned char* p){unsigned v=0;for(int i=0;i<4;i++)v|=unsigned(p[i])<<(8*i);return v;}
void require(bool b,const char* s){if(!b)throw std::runtime_error(s);}
Bytes unpack(const Bytes& f,size_t header,U n,U size,size_t hash_at){
 require(f.size()==header+n && size<=128*1024*1024,"frame sizes");
 Bytes b(size);z_stream z{};z.next_in=const_cast<unsigned char*>(f.data()+header);z.avail_in=n;z.next_out=b.data();z.avail_out=size;
 require(inflateInit(&z)==Z_OK,"inflate init");int rc=inflate(&z,Z_FINISH);bool ok=rc==Z_STREAM_END && z.total_in==n && z.total_out==size;inflateEnd(&z);require(ok,"inflate framing");
 std::array<unsigned char,32> h{};SHA256(b.data(),b.size(),h.data());require(std::equal(h.begin(),h.end(),f.begin()+hash_at),"frame SHA256");return b;
}
Bytes read(std::ifstream& f,const std::string& dir,const std::string& name,std::string& current,U off,U size){
 if(name!=current){if(f.is_open())f.close();f.open(dir+"/"+name,std::ios::binary);current=name;require(bool(f),"pack open");}
 f.seekg(off);Bytes b(size);require(bool(f.read(reinterpret_cast<char*>(b.data()),size)),"pack read");return b;
}
U var(const Bytes& b,size_t& at){U x=0;for(unsigned j=0;j<10;j++){require(at<b.size(),"varint end");unsigned c=b[at++];require(j<9 || c<=1,"varint overflow");x|=U(c&127)<<(7*j);if(c<128){require(!j || c!=0,"varint canonical");return x;}}throw std::runtime_error("varint length");}
U power(U a,U e,U p){U y=1;while(e){if(e&1)y=y*a%p;a=a*a%p;e>>=1;}return y;}
const int small_primes[]={3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127};
struct Acc {
 U n=0,edge=0,tie=0,eqML=0;long double m=0,m2=0,l=0,k=0,z=0,z2=0,mz=0,a=0,a2=0,ma=0;
 void add(U M,U L,U K,long double Z,U A,U boundary=0){n++;m+=M;m2+=M*M;l+=L;k+=K;z+=Z;z2+=Z*Z;mz+=M*Z;a+=A;a2+=static_cast<long double>(A)*A;ma+=static_cast<long double>(M)*A;edge+=K>M+1;tie+=K==M+1 && K==boundary;eqML+=M==L;}
 void out(std::ostream& o)const{o<<"{\"n\":"<<n<<",\"edge\":"<<edge<<",\"tie\":"<<tie<<",\"eqML\":"<<eqML<<",\"m\":"<<m<<",\"m2\":"<<m2<<",\"l\":"<<l<<",\"k\":"<<k<<",\"z\":"<<z<<",\"z2\":"<<z2<<",\"mz\":"<<mz<<",\"a\":"<<a<<",\"a2\":"<<a2<<",\"ma\":"<<ma<<"}";}
};
struct Cross {U n=0,eq=0;long double x=0,y=0,xx=0,yy=0,xy=0;void add(U a,U b){n++;eq+=a==b;x+=a;y+=b;xx+=a*a;yy+=b*b;xy+=a*b;}void out(std::ostream&o)const{o<<"{\"n\":"<<n<<",\"eq\":"<<eq<<",\"x\":"<<x<<",\"y\":"<<y<<",\"xx\":"<<xx<<",\"yy\":"<<yy<<",\"xy\":"<<xy<<"}";}};
struct Gap {U n=0,first=0,first_rank=0,last=0,last_rank=0,max_rank=0,rank_left=0,rank_right=0,max_p=0,p_left=0,p_right=0;
 void add(U p,U rank){if(n==0){first=p;first_rank=rank;}else{U g=rank-last_rank-1;if(g>max_rank){max_rank=g;rank_left=last;rank_right=p;}if(p-last>max_p){max_p=p-last;p_left=last;p_right=p;}}last=p;last_rank=rank;n++;}
 void out(std::ostream&o,U eligible)const{o<<"{\"count\":"<<n<<",\"first\":"<<first<<",\"last\":"<<last<<",\"max_rank_gap\":"<<max_rank<<",\"rank_left\":"<<rank_left<<",\"rank_right\":"<<rank_right<<",\"max_integer_gap\":"<<max_p<<",\"integer_left\":"<<p_left<<",\"integer_right\":"<<p_right<<",\"left_censored\":"<<(n?first_rank-1:eligible)<<",\"right_censored\":"<<(n?eligible-last_rank:eligible)<<"}";}
};
struct Target {int r,k;U eligible=0;Gap a,b;U bands[16][2]{};};
struct Record {U p=0,m=0,l=0,k=0,w=0,c=0,a=0;int r=0,s=0;void out(std::ostream&o)const{o<<"{\"r\":"<<r<<",\"p\":"<<p<<",\"M\":"<<m<<",\"L\":"<<l<<",\"K\":"<<k<<",\"witness\":"<<w<<",\"sigma\":"<<s<<",\"C_joined\":"<<c<<",\"A4\":"<<a<<"}";}};
void keep(std::vector<Record>& v,Record q,int which){v.push_back(q);std::sort(v.begin(),v.end(),[which](const Record&a,const Record&b){U x=which==0?a.m:which==1?a.l:a.k,y=which==0?b.m:which==1?b.l:b.k;return x!=y?x>y:a.p>b.p;});if(v.size()>12)v.resize(12);}
Acc stats[8][17],classes[8][4],traceclass[24],tracebands[16],traceM[128];
Cross cross[8][8][17];
U hist[8][3][256]{},bands[8][16][128]{},joint[8][128][128]{},by_sigma[8][2][128]{},bestM[8][128]{},bestL[8][128]{};
U tailN[8]{},tailCenter[8]{},tailPrefix[8]{},tailWitness[8]{};
int cuts[8]={40,28,23,20,19,18,17,16};
U zhist[80]{},totalP=0,sumP=0,sumA=0,traceN=0,blocks=0,bytes=0;
long double sumZ=0,sumZ2=0,sumRatio=0,maxabsZ=0,maxRelative=0;I sumT=0;U maxZP=0,maxRelP=0;
std::vector<Record> records[8][3];
int main(int argc,char**argv){try{
 require(argc==5,"usage: stream scan_dir A4_dir manifest output");
 auto began=std::chrono::steady_clock::now();std::ifstream mf(argv[3]),sf,af;std::string sn,an,cs,ca;U id,lo,hi,count,so,ss,ao,as;
 std::vector<Target> targets;for(int k=25;k<=28;k++)targets.push_back(Target{2,k});for(int k=17;k<=25;k++)targets.push_back(Target{3,k});
 U expected_lo=1;double maxblock=0;
 while(mf>>id>>lo>>hi>>count>>sn>>so>>ss>>an>>ao>>as){auto bt=std::chrono::steady_clock::now();require(id==blocks && lo==expected_lo && hi>lo,"coverage");
  Bytes f=read(sf,argv[1],sn,cs,so,ss),g=read(af,argv[2],an,ca,ao,as);bytes+=ss+as;
  require(f.size()>=80 && std::string(reinterpret_cast<char*>(f.data()),8)==std::string("VDWSM1\0\0",8),"scan magic");
  require(u64(f.data()+24)==lo && u64(f.data()+32)==hi && u32(f.data()+40)==count && u32(f.data()+44)==1,"scan identity");
  Bytes raw=unpack(f,80,u64(f.data()+8),u64(f.data()+16),48);
  require(g.size()>=136 && std::string(reinterpret_cast<char*>(g.data()),8)==std::string("VDWA4V1\0",8),"A4 magic");
  require(u64(g.data()+8)==lo && u64(g.data()+16)==hi && u32(g.data()+24)==count && u32(g.data()+36)==0 && u32(g.data()+32)==count*4,"A4 identity");
  Bytes av=unpack(g,136,u32(g.data()+28),u32(g.data()+32),104);
  size_t at=0;U p=lo-1;
  for(U row=0;row<count;row++){
   U d=var(raw,at);require(d>0,"prime delta");p+=d;require(p>=lo && p<hi && at<raw.size(),"prime range");
   unsigned mask=raw[at++],expected_mask=0;for(int r=2;r<=9;r++)if((p-1)%r==0)expected_mask|=1<<(r-2);require(mask==expected_mask,"mask");
   U A=u32(av.data()+4*row);if(p==2)require(A==0,"A4 placeholder");else sumA+=A;
   totalP++;sumP+=p;U m[8]{},l[8]{},k[8]{},w[8]{},c[8]{};int sigma[8]{};int band=p>>28;require(band<16,"band");
   long double Z=0;
   if(p>3){int c2=(p%8==1 || p%8==7)?1:-1,c3=(p%12==1 || p%12==11)?1:-1,cm=p%4==1?1:-1;int B=4*c2+2*c2*c3+2*c3+2*cm+2*cm*c2;
    I T=I(p)-11-B-8*I(A);U radius=std::sqrt(static_cast<long double>(4*p));while((radius+1)*(radius+1)<=4*p)radius++;while(radius*radius>4*p)radius--;require(std::abs(T)<=I(radius),"Hasse bound");
    Z=T/(2*std::sqrt(static_cast<long double>(p)));traceN++;sumT+=T;sumZ+=Z;sumZ2+=Z*Z;sumRatio+=static_cast<long double>(A)/p;
    int zi=std::clamp(int((Z+1)*40),0,79);zhist[zi]++;if(std::abs(Z)>maxabsZ){maxabsZ=std::abs(Z);maxZP=p;}
    if(p>=1000000){long double re=std::abs(static_cast<long double>(8)*A/p-1);if(re>maxRelative){maxRelative=re;maxRelP=p;}}
   }
   for(int i=0;i<8;i++)if(mask&(1<<i)){
    int r=i+2;m[i]=var(raw,at);l[i]=var(raw,at);w[i]=var(raw,at);sigma[i]=((p-1)/2)%r==0?1:-1;
    require(1<=l[i] && l[i]<=m[i] && m[i]<=(p-1)/r && w[i]>=1 && w[i]<=(p-1)/2,"profile bounds");
    k[i]=std::max(m[i]+1,sigma[i]==1?2*l[i]+2:l[i]+2);require(m[i]<128 && l[i]<128 && k[i]<256,"analysis histogram capacity");
    hist[i][0][m[i]]++;hist[i][1][l[i]]++;hist[i][2][k[i]]++;bands[i][band][m[i]]++;joint[i][m[i]][l[i]]++;by_sigma[i][sigma[i]==1][m[i]]++;
    bestM[i][m[i]]=p;bestL[i][l[i]]=p;stats[i][16].add(m[i],l[i],k[i],Z,A,sigma[i]==1?2*l[i]+2:l[i]+2);stats[i][band].add(m[i],l[i],k[i],Z,A,sigma[i]==1?2*l[i]+2:l[i]+2);
    classes[i][(sigma[i]==1)*2+(l[i]>1)].add(m[i],l[i],k[i],Z,A,sigma[i]==1?2*l[i]+2:l[i]+2);
    if(m[i]>=U(cuts[i])){tailN[i]++;tailPrefix[i]+=l[i]==m[i];if(sigma[i]==1){
      if(l[i]>1)c[i]=l[i];else for(int q:small_primes){require(U(q)<=m[i]+1,"central upper bound");if(power(q,(p-1)/r,p)!=1){c[i]=q-1;break;}}
      require(c[i]>0 && c[i]<=m[i],"central length");if(c[i]==m[i]){tailCenter[i]++;tailWitness[i]+=w[i]==(p-c[i]+1)/2;}
    }}
    Record rec{p,m[i],l[i],k[i],w[i],c[i],A,r,sigma[i]};
    if(records[i][0].size()<12 || m[i]>=records[i][0].back().m)keep(records[i][0],rec,0);
    if(records[i][1].size()<12 || l[i]>=records[i][1].back().l)keep(records[i][1],rec,1);
    if(records[i][2].size()<12 || k[i]>=records[i][2].back().k)keep(records[i][2],rec,2);
   }
   if(p>3){traceclass[p%24].add(m[0],l[0],k[0],Z,A);tracebands[band].add(m[0],l[0],k[0],Z,A);traceM[m[0]].add(m[0],l[0],k[0],Z,A);}
   for(int i=0;i<8;i++)if(m[i])for(int j=i+1;j<8;j++)if(m[j]){
    if((j+2)%(i+2)==0)require(m[j]<=m[i] && l[j]<=l[i] && k[j]<=k[i],"refinement inequality");
    cross[i][j][16].add(m[i],m[j]);cross[i][j][band].add(m[i],m[j]);
   }
   for(auto& t:targets)if(m[t.r-2] && p>U(t.k)){int i=t.r-2;t.eligible++;if(m[i]<U(t.k)){t.a.add(p,t.eligible);t.bands[band][0]++;}if(k[i]<=U(t.k)){t.b.add(p,t.eligible);t.bands[band][1]++;}}
  }
  require(at==raw.size(),"trailing bytes");expected_lo=hi;blocks++;double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-bt).count();maxblock=std::max(maxblock,sec);require(sec<30,"analysis batch deadline");
  if(blocks%2048==0)std::cerr<<"blocks="<<blocks<<" primes="<<totalP<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<'\n';
 }
 std::ofstream o(argv[4]);o<<std::setprecision(19)<<"{\"blocks\":"<<blocks<<",\"next_lo\":"<<expected_lo<<",\"primes\":"<<totalP<<",\"sum_p\":"<<sumP<<",\"sum_A4\":"<<sumA<<",\"bytes_read\":"<<bytes<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<",\"max_block_seconds\":"<<maxblock;
 o<<",\"hist\":{";for(int i=0;i<8;i++){if(i)o<<',';o<<'"'<<i+2<<"\":[";for(int j=0;j<3;j++){if(j)o<<',';o<<'{';bool comma=false;for(int v=0;v<256;v++)if(hist[i][j][v]){if(comma)o<<',';comma=true;o<<'"'<<v<<"\":"<<hist[i][j][v];}o<<'}';}o<<']';}o<<'}';
 o<<",\"stats\":{";for(int i=0;i<8;i++){if(i)o<<',';o<<'"'<<i+2<<"\":[";for(int b=0;b<17;b++){if(b)o<<',';stats[i][b].out(o);}o<<']';}o<<'}';
 o<<",\"classes\":{";for(int i=0;i<8;i++){if(i)o<<',';o<<'"'<<i+2<<"\":[";for(int b=0;b<4;b++){if(b)o<<',';classes[i][b].out(o);}o<<']';}o<<'}';
 o<<",\"band_hist_M\":{";for(int i=0;i<8;i++){if(i)o<<',';o<<'"'<<i+2<<"\":[";for(int b=0;b<16;b++){if(b)o<<',';o<<'{';bool comma=false;for(int v=0;v<128;v++)if(bands[i][b][v]){if(comma)o<<',';comma=true;o<<'"'<<v<<"\":"<<bands[i][b][v];}o<<'}';}o<<']';}o<<'}';
 o<<",\"M_by_sigma\":{";for(int i=0;i<8;i++){if(i)o<<',';o<<'"'<<i+2<<"\":[";for(int b=0;b<2;b++){if(b)o<<',';o<<'{';bool comma=false;for(int v=0;v<128;v++)if(by_sigma[i][b][v]){if(comma)o<<',';comma=true;o<<'"'<<v<<"\":"<<by_sigma[i][b][v];}o<<'}';}o<<']';}o<<'}';
 o<<",\"joint_ML\":{";for(int i=0;i<8;i++){if(i)o<<',';o<<'"'<<i+2<<"\":[";bool comma=false;for(int m=0;m<128;m++)for(int l=0;l<128;l++)if(joint[i][m][l]){if(comma)o<<',';comma=true;o<<'['<<m<<','<<l<<','<<joint[i][m][l]<<']';}o<<']';}o<<'}';
 o<<",\"best_M\":{";for(int i=0;i<8;i++){if(i)o<<',';o<<'"'<<i+2<<"\":{";bool comma=false;for(int v=0;v<128;v++)if(bestM[i][v]){if(comma)o<<',';comma=true;o<<'"'<<v<<"\":"<<bestM[i][v];}o<<'}';}o<<'}';
 o<<",\"best_L\":{";for(int i=0;i<8;i++){if(i)o<<',';o<<'"'<<i+2<<"\":{";bool comma=false;for(int v=0;v<128;v++)if(bestL[i][v]){if(comma)o<<',';comma=true;o<<'"'<<v<<"\":"<<bestL[i][v];}o<<'}';}o<<'}';
 o<<",\"cross\":{";bool comma=false;for(int i=0;i<8;i++)for(int j=i+1;j<8;j++){if(comma)o<<',';comma=true;o<<'"'<<i+2<<','<<j+2<<"\":[";for(int b=0;b<17;b++){if(b)o<<',';cross[i][j][b].out(o);}o<<']';}o<<'}';
 o<<",\"tails\":{";for(int i=0;i<8;i++){if(i)o<<',';o<<'"'<<i+2<<"\":{\"cut\":"<<cuts[i]<<",\"n\":"<<tailN[i]<<",\"central_attains_M\":"<<tailCenter[i]<<",\"prefix_attains_M\":"<<tailPrefix[i]<<",\"witness_at_center\":"<<tailWitness[i]<<'}';}o<<'}';
 o<<",\"records\":{";for(int i=0;i<8;i++){if(i)o<<',';o<<'"'<<i+2<<"\":[";for(int j=0;j<3;j++){if(j)o<<',';o<<'[';for(size_t a=0;a<records[i][j].size();a++){if(a)o<<',';records[i][j][a].out(o);}o<<']';}o<<']';}o<<'}';
 o<<",\"targets\":{";comma=false;for(const auto&t:targets){if(comma)o<<',';comma=true;o<<'"'<<t.r<<','<<t.k<<"\":{\"eligible\":"<<t.eligible<<",\"A\":";t.a.out(o,t.eligible);o<<",\"AB\":";t.b.out(o,t.eligible);o<<",\"bands\":[";for(int b=0;b<16;b++){if(b)o<<',';o<<'['<<t.bands[b][0]<<','<<t.bands[b][1]<<']';}o<<"]}";}o<<'}';
 o<<",\"A4\":{\"trace_n\":"<<traceN<<",\"sum_t\":"<<sumT<<",\"sum_z\":"<<sumZ<<",\"sum_z2\":"<<sumZ2<<",\"sum_A4_over_p\":"<<sumRatio<<",\"max_abs_z\":"<<maxabsZ<<",\"max_abs_z_p\":"<<maxZP<<",\"max_relative_error_p_ge_1e6\":"<<maxRelative<<",\"max_relative_error_p\":"<<maxRelP<<",\"trace_hist\":[";
 for(int b=0;b<80;b++){if(b)o<<',';o<<zhist[b];}o<<"],\"mod24\":[";for(int b=0;b<24;b++){if(b)o<<',';traceclass[b].out(o);}o<<"],\"bands\":[";for(int b=0;b<16;b++){if(b)o<<',';tracebands[b].out(o);}o<<"],\"by_M2\":[";for(int b=0;b<128;b++){if(b)o<<',';traceM[b].out(o);}o<<"]}}\n";require(bool(o),"output write");
 std::cerr<<"COMPLETE blocks="<<blocks<<" primes="<<totalP<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
