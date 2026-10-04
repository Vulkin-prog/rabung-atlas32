#define main unused_archive_reader_main
#include "../analyse_atlas_32_20261004/stream.cpp"
#undef main
static const int pairs[5][2]={{3,9},{4,8},{3,6},{2,4},{2,3}};
static const int deep_cuts[5][2]={{26,20},{24,20},{30,20},{46,26},{46,30}};
static const int needed_cut[10]={0,0,46,26,24,128,20,128,20,20};
static U jm[5][16][16][128][128]{};
static U jk[5][16][256][256]{};
static U jt[5][16][16][3][3]{};
int main(int argc,char**argv){try{
 require(argc==4,"usage scan_dir frames output");
 auto began=std::chrono::steady_clock::now();std::ifstream mf(argv[2]),sf;
 std::string name,an,current;U id,lo,hi,count,off,size,ao,as,nb=0,np=0,next=1;
 std::vector<Record> residuals;U mh[8][128]{},profile_count[8]{},tail_count[2]{},tail_union[2]{};
 double maximum_block=0;
 while(mf>>id>>lo>>hi>>count>>name>>off>>size>>an>>ao>>as){
  auto bt=std::chrono::steady_clock::now();require(id==nb&&lo==next,"coverage");
  Bytes f=read(sf,argv[1],name,current,off,size);
  require(f.size()>=80&&std::string(reinterpret_cast<char*>(f.data()),8)==std::string("VDWSM1\0\0",8),"frame magic");
  require(u64(f.data()+24)==lo&&u64(f.data()+32)==hi&&u32(f.data()+40)==count&&u32(f.data()+44)==1,"frame identity");
  Bytes raw=unpack(f,80,u64(f.data()+8),u64(f.data()+16),48);U p=lo-1;size_t at=0;
  for(U row=0;row<count;row++){
   U delta=var(raw,at);require(delta>0,"positive prime delta");p+=delta;require(p>=lo&&p<hi,"prime range");unsigned mask=raw.at(at++),expected=0;
   U m[10]{},l[10]{},k[10]{},w[10]{},center[10]{};bool sig[10]{};
   for(int r=2;r<=9;r++)if((p-1)%r==0)expected|=1u<<(r-2);require(mask==expected,"mask");
   for(int r=2;r<=9;r++)if(mask&(1u<<(r-2))){
    m[r]=var(raw,at);l[r]=var(raw,at);w[r]=var(raw,at);sig[r]=((p-1)/2)%r==0;
    k[r]=std::max(m[r]+1,sig[r]?2*l[r]+2:l[r]+2);
    require(m[r]<128&&k[r]<256,"hist bounds");mh[r-2][m[r]]++;profile_count[r-2]++;
    if(sig[r]&&m[r]>=needed_cut[r]){
     if(l[r]>1)center[r]=l[r];else for(int q:small_primes)if(power(q,(p-1)/r,p)!=1){center[r]=q-1;break;}
     require(center[r]>0&&center[r]<=m[r],"needed center check");
    }
    if((r==8||r==9)&&m[r]>=20){
     tail_count[r-8]++;U c=center[r];
     if(c==m[r]||l[r]==m[r])tail_union[r-8]++;
     else residuals.push_back(Record{p,m[r],l[r],k[r],w[r],c,0,r,sig[r]?1:-1});
    }
   }
   int band=p>>28;
   for(int q=0;q<5;q++){
    int a=pairs[q][0],b=pairs[q][1];if(!m[a]||!m[b])continue;
    if(b%a==0)require(m[b]<=m[a]&&l[b]<=l[a]&&k[b]<=k[a],"refinement");
    unsigned cls=(l[a]>1)+2*(l[b]>1)+4*sig[a]+8*sig[b];
    jm[q][band][cls][m[a]][m[b]]++;jk[q][band][k[a]][k[b]]++;
    int ta=m[a]<deep_cuts[q][0]?0:(center[a]==m[a]||l[a]==m[a])?1:2;
    int tb=m[b]<deep_cuts[q][1]?0:(center[b]==m[b]||l[b]==m[b])?1:2;
    jt[q][band][cls][ta][tb]++;
   }
   np++;
  }
  require(at==raw.size(),"trailing bytes");nb++;next=hi;
  double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-bt).count();
  maximum_block=std::max(maximum_block,sec);require(sec<30,"batch deadline");
  if(nb%8192==0)std::cerr<<"blocks="<<nb<<" primes="<<np<<'\n';
 }
 require(nb==49901&&np==203280221&&next==(1ull<<32),"full coverage");
 require(tail_count[0]==35&&tail_count[1]==21&&tail_union[0]==12&&tail_union[1]==7&&residuals.size()==37,"residual reconciliation");
 std::ofstream o(argv[3]);o<<std::setprecision(15)<<"{\"status\":\"PASS\",\"blocks\":"<<nb<<",\"primes\":"<<np
 <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<",\"max_block_seconds\":"<<maximum_block
 <<",\"class_bits\":[\"L_coarse>1\",\"L_fine>1\",\"sigma_coarse=1\",\"sigma_fine=1\"],\"residuals\":[";
 for(size_t i=0;i<residuals.size();i++){if(i)o<<',';auto z=residuals[i];o<<"{\"r\":"<<z.r<<",\"p\":"<<z.p<<",\"M\":"<<z.m<<",\"L\":"<<z.l<<",\"K\":"<<z.k<<",\"witness\":"<<z.w<<",\"sigma\":"<<z.s<<",\"C_joined\":"<<z.c<<"}";}o<<"],\"M_histograms\":{";
 for(int r=2;r<=9;r++){if(r>2)o<<',';o<<'"'<<r<<"\":{";bool comma=false;for(int m=0;m<128;m++)if(mh[r-2][m]){if(comma)o<<',';comma=true;o<<'"'<<m<<"\":"<<mh[r-2][m];}o<<'}';}o<<"},\"pairs\":{";
 for(int q=0;q<5;q++){if(q)o<<',';o<<'"'<<pairs[q][0]<<','<<pairs[q][1]<<"\":{\"M_by_band_class\":[";bool comma=false;
  for(int b=0;b<16;b++)for(int c=0;c<16;c++)for(int a=0;a<128;a++)for(int d=0;d<128;d++)if(jm[q][b][c][a][d]){
   if(comma)o<<',';comma=true;o<<'['<<b<<','<<c<<','<<a<<','<<d<<','<<jm[q][b][c][a][d]<<']';}
  o<<"],\"K_by_band\":[";comma=false;for(int b=0;b<16;b++)for(int a=0;a<256;a++)for(int d=0;d<256;d++)if(jk[q][b][a][d]){
   if(comma)o<<',';comma=true;o<<'['<<b<<','<<a<<','<<d<<','<<jk[q][b][a][d]<<']';}
  o<<"],\"deep_cuts\":["<<deep_cuts[q][0]<<','<<deep_cuts[q][1]<<"],\"tail_status_labels\":[\"below_cut\",\"center_or_prefix_attains_M\",\"residual\"],\"tail_by_band_class\":[";comma=false;
  for(int b=0;b<16;b++)for(int c=0;c<16;c++)for(int a=0;a<3;a++)for(int d=0;d<3;d++)if(jt[q][b][c][a][d]){
   if(comma)o<<',';comma=true;o<<'['<<b<<','<<c<<','<<a<<','<<d<<','<<jt[q][b][c][a][d]<<']';}o<<"]}";
 }
 o<<"}}\n";require(bool(o),"output write");return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
