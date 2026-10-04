#define main unused_archive_reader_main
#include "../analyse_atlas_32_20261004/stream.cpp"
#undef main
#include <numeric>
// Counts the exact phase of the canonical archived maximum, all r8/r9 profiles.
// Classes: L>1; sigma=1. Deep classification only is needed above M>=13.
static U phase_hist[2][128][4][2]{};
static U deep_hist[2][128][4][2][3]{};
static U quad_hist[2][128][4][2][2]{};
int main(int argc,char**argv){try{
 require(argc==4,"usage scan_dir frames output");
 auto began=std::chrono::steady_clock::now();std::ifstream mf(argv[2]),sf;
 std::string name,an,current;U id,lo,hi,count,off,size,ao,as,nb=0,np=0,next=1;
 double max_block=0;
 while(mf>>id>>lo>>hi>>count>>name>>off>>size>>an>>ao>>as){
  auto bt=std::chrono::steady_clock::now();require(id==nb&&lo==next,"coverage"); if(id!=0&&id!=1&&id!=1024&&id!=4096&&id!=8192&&id!=16384&&id!=32768&&id!=49900){nb++;next=hi;continue;}
  Bytes f=read(sf,argv[1],name,current,off,size);
  require(f.size()>=80&&std::string(reinterpret_cast<char*>(f.data()),8)==std::string("VDWSM1\0\0",8),"magic");
  require(u64(f.data()+24)==lo&&u64(f.data()+32)==hi&&u32(f.data()+40)==count&&u32(f.data()+44)==1,"identity");
  Bytes raw=unpack(f,80,u64(f.data()+8),u64(f.data()+16),48);U p=lo-1;size_t at=0;
  for(U row=0;row<count;row++){
   U delta=var(raw,at);require(delta>0,"delta");p+=delta;require(p>=lo&&p<hi,"prime range");
   unsigned mask=raw.at(at++),expected=0;for(int r=2;r<=9;r++)if((p-1)%r==0)expected|=1u<<(r-2);require(mask==expected,"mask");
   for(int r=2;r<=9;r++)if(mask&(1u<<(r-2))){
    U m=var(raw,at),l=var(raw,at),w=var(raw,at);if(r<8)continue;
    require(m<128&&w>=1&&w+m<=p,"profile");
    bool sig=((p-1)/2)%r==0;int cls=(l>1)+2*sig,neutral=power(w,(p-1)/r,p)==1;
    phase_hist[r-8][m][cls][neutral]++;
    if(m>=13){
     U c=0;if(sig){if(l>1)c=l;else for(int q:small_primes)if(power(q,(p-1)/r,p)!=1){c=q-1;break;}require(c>0&&c<=m,"center");}
     int category=c==m||l==m?0:2;
     if(category==2){
      for(U q=2;q<=128&&category==2;q++){
       U a0=(q*w+p-1)/p,a1=q*(w+m-1)/p;
       for(U a=a0;a<=a1;a++)if(std::gcd(a,q)==1){category=1;break;}
      }
     }
     deep_hist[r-8][m][cls][neutral][category]++;
     if(category==2){bool hit=false;for(I A=-64;A<=64&&!hit;A++){I t=(2*(I)w+A)%(I)p;if(t<0)t+=p;U square=U(t)*U(t)%p;I D=(I)square;if(D>(I)p/2)D-=p;if(D==0||D< -256||D>256||(A*A-D)%4)continue;bool rational=false;for(I k=1;k<=16;k++)if(D==k*k)rational=true;if(!rational)hit=true;}quad_hist[r-8][m][cls][neutral][hit]++;}
    }
   }
   np++;
  }
  require(at==raw.size(),"trailing bytes");nb++;next=hi;
  double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-bt).count();max_block=std::max(max_block,sec);require(sec<30,"batch deadline");
  if(nb%8192==0)std::cerr<<"blocks="<<nb<<" primes="<<np<<'\n';
 }
 require(nb==49901&&next==(1ull<<32),"complete coverage");
 std::ofstream o(argv[3]);o<<std::setprecision(15)<<"{\"status\":\"PASS\",\"blocks\":"<<nb<<",\"primes\":"<<np<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<",\"max_block_seconds\":"<<max_block
 <<",\"class_bits\":[\"L>1\",\"sigma=1\"],\"deep_category_labels\":[\"center_or_prefix_attains_M\",\"rational_cut_q_le_128_inside_canonical_maximum\",\"residual\"],\"counts\":[";
 bool comma=false;for(int r=8;r<=9;r++)for(int m=0;m<128;m++)for(int c=0;c<4;c++)for(int n=0;n<2;n++)if(phase_hist[r-8][m][c][n]){
  if(comma)o<<',';comma=true;o<<'['<<r<<','<<m<<','<<c<<','<<n<<','<<phase_hist[r-8][m][c][n]<<']';}
 o<<"],\"deep_counts\":[";comma=false;for(int r=8;r<=9;r++)for(int m=13;m<128;m++)for(int c=0;c<4;c++)for(int n=0;n<2;n++)for(int k=0;k<3;k++)if(deep_hist[r-8][m][c][n][k]){
  if(comma)o<<',';comma=true;o<<'['<<r<<','<<m<<','<<c<<','<<n<<','<<k<<','<<deep_hist[r-8][m][c][n][k]<<']';}o<<"],\"quadratic_detection\":{\"A_bound\":64,\"D_bound\":256,\"nonzero_nonsquare_D\":true},\"quadratic_counts\":[";comma=false;for(int r=8;r<=9;r++)for(int m=13;m<128;m++)for(int c=0;c<4;c++)for(int n=0;n<2;n++)for(int k=0;k<2;k++)if(quad_hist[r-8][m][c][n][k]){if(comma)o<<',';comma=true;o<<'['<<r<<','<<m<<','<<c<<','<<n<<','<<k<<','<<quad_hist[r-8][m][c][n][k]<<']';}o<<"]}\n";require(bool(o),"output");return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
