#define VDW_STREAM_NO_MAIN
#include "../reference/stream_reference.cpp"
#include "positive_filter.hpp"
#include "wheel_filter.hpp"
#include "sparse_words.hpp"
struct Positive {
 const Table&T;u32 stride,np,flip;std::vector<u32> data;
 explicit Positive(const Table&t):T(t),stride((t.B+32)/32+3),np(t.planes.size()),flip(0){
  u32 two=1;while(T.r%(two*2)==0)two*=2;if(T.sigma)flip=two/2;
  data.assign(np*stride,0);
  for(u32 x=1;x<=T.B;++x){u32 e=T.c[x]%two+two*(T.c[x]%(T.r/two));
   for(u32 j=0;j<np;++j)if(e&(1u<<j))data[j*stride+1+(x>>5)]|=1u<<(x&31);}
 }
 vdw_next::View view()const{return {data.data(),stride,T.B,np,flip};}
};
struct Wheel {
 const Table&T;u32 g,np,stride,flip;std::vector<u32> residues,data;
 Wheel(const Table&t,u32 gg):T(t),g(gg),np(t.planes.size()),stride(((t.B/g)+32)/32+3),flip(0){
  require(g==2||g==3||g==6,"wheel modulus");u32 two=1;while(T.r%(two*2)==0)two*=2;if(T.sigma)flip=two/2;
  for(u32 c=1;c<g;++c)if(std::gcd(c,g)==1)residues.push_back(c);
  data.assign(residues.size()*np*stride,0);
  for(u32 row=0;row<residues.size();++row)for(u32 b=residues[row],n=0;b<=T.B;b+=g,++n){
   u32 e=T.c[b]%two+two*(T.c[b]%(T.r/two));for(u32 k=0;k<np;++k)if(e&(1u<<k))data[(row*np+k)*stride+1+(n>>5)]|=1u<<(n&31);}
 }
 vdw_next::ResidueView view(u32 row)const{
  u32 c=residues.at(row),mr=0;while(residues.at(mr)!=g-c)++mr;
  return {data.data()+row*np*stride,data.data()+mr*np*stride,stride,np,flip,(T.B-c)/g,(T.B-(g-c))/g};
 }
};
struct NewResult {
 StreamResult result;vdw_next::Counters filterstats;
 u64 wheel_words[7]={},wheel_plane_bytes[7]={};
 u64 sparse_words=0,nonempty_warps=0,output_bits_before_gcd=0;
 std::vector<vdw_next::SparseWord> descriptors;
 u64 words=0,relaxed_bits=0,primitive=0,survivor_reps=0,orientation_total=0;
 u32 seed=0,seed_witness=0,target=0;
 std::vector<std::pair<u32,u32>> representatives;std::vector<std::tuple<u32,u32,u32>> potential;
};
static NewResult improved(const Table&T,const Positive&P,bool pruning=true,bool audit=false,bool instrument=true,bool wheel=false){
 NewResult N;auto &R=N.result;R.multiplicity_bound=mult_bound(T);R.upper=select_t(T.p,T.r).upper;
 u32 run=0;int prev=-1;
 for(u32 x=1;x<=T.B;++x){int c=T.c[x];run=c==prev?run+1:1;prev=c;
  if(run>R.M){R.M=run;R.witness=x-run+1;}}
 N.seed=R.M;N.seed_witness=R.witness;N.target=std::min(256u,std::max(T.L,R.M+1));
 while(R.prefix<T.B&&T.c[R.prefix+1]==0)++R.prefix;
 // A table-constant prefix must be completed even if a future filter excludes it.
 if(R.prefix==T.B){while((u64)R.prefix+1<T.p&&small_oracle(T,R.prefix+1,R.st)==0)++R.prefix;
  if(R.prefix>R.M){R.M=R.prefix;R.witness=1;}}
 std::vector<u32> inv(T.D+1);inv[1]=1;for(u32 d=2;d<=T.D;++d)inv[d]=T.p-(u64)(T.p/d)*inv[T.p%d]%T.p;
 auto visit=[&](u32 x,i64 a,u32 d){++R.st.oriented;
  require(x>=1&&(u64)x+T.L<=T.p,"new linear lift");
  const int col=(T.colour(a)+T.r-T.c[d])%T.r;
  auto at=[&](i64 z,u32 xx){return z&&z>=-(i64)T.B&&z<=T.B?(T.colour(z)+T.r-T.c[d])%T.r:small_oracle(T,xx,R.st);};
  if(x>1&&at(a-(i64)d,x-1)==col)return;
  ++R.st.head_visits;u32 len=T.L;
  while((u64)x+len<T.p){int c=at(a+(i64)len*d,x+len);++R.st.extensions;if(c!=col)break;++len;}
  if(len>R.M||(len==R.M&&x<R.witness)){R.M=len;R.witness=x;}
  if(x==1)R.prefix=len;
 };
 Wheel W2(T,2),W3(T,3),W6(T,6);
 N.wheel_plane_bytes[1]=P.data.size()*4;N.wheel_plane_bytes[2]=W2.data.size()*4;N.wheel_plane_bytes[3]=W3.data.size()*4;N.wheel_plane_bytes[6]=W6.data.size()*4;
 auto traverse=[&](auto view,u32 d,u32 g,u32 c){
  const u32 minimum=((T.L-1)*d+1)/2,first=minimum>c?(minimum-c+g-1)/g:0,maximum=(T.B-c)/g,step=d/g;
  if(first>maximum)return;
  u32 last_nonempty_warp=UINT32_MAX;
  for(u32 w=first>>5;w<=maximum/32;++w){++N.words;++N.wheel_words[g];
   u32 bits=vdw_next::filter_range(view,step,T.L,w,first,maximum,instrument?&N.filterstats:nullptr);N.relaxed_bits+=__builtin_popcount(bits);
   auto keep=pruning?vdw_next::prefilter_row(view,step,T.L,N.target,w,bits,instrument?&N.filterstats:nullptr):vdw_next::Oriented{bits,bits};
   u32 rest=keep.left|keep.right;
   if(rest){++N.sparse_words;N.output_bits_before_gcd+=__builtin_popcount(rest);
    u32 warp=(w-(first>>5))/32;if(warp!=last_nonempty_warp){++N.nonempty_warps;last_nonempty_warp=warp;}
    if(audit){u32 tag=g==1?0:(g==2?1:(g==3?2:3));N.descriptors.push_back({d|(tag<<30),g*32*w+c,keep.left,keep.right});}
   }
   if(audit){for(u32 z=bits;z;z&=z-1){u32 b=g*(32*w+__builtin_ctz(z))+c;if(vdw_next::admissible(b,d,T.L))N.representatives.emplace_back(d,b);}}
   while(rest){u32 bit=__builtin_ctz(rest);rest&=rest-1;u32 b=g*(32*w+bit)+c;if(!vdw_next::admissible(b,d,T.L))continue;
    ++N.survivor_reps;u32 flags=((keep.left>>bit)&1)|(((keep.right>>bit)&1)<<1);N.orientation_total+=__builtin_popcount(flags);
    if(audit)N.potential.emplace_back(d,b,flags);
    u32 x=(T.p-(u64)b*inv[d]%T.p)%T.p,xr=T.p-x-(T.L-1);
    if(flags&1)visit(x,-(i64)b,d);
    if((flags&2)&&(!(flags&1)||xr!=x))visit(xr,(i64)b-(T.L-1)*(i64)d,d);
   }
  }
 };
 for(u32 d=1;d<=T.D;++d){u32 g=wheel?std::gcd(d,6u):1;
  if(g==1)traverse(P.view(),d,1,0);
  else{const Wheel&W=g==2?W2:(g==3?W3:W6);for(u32 row=0;row<W.residues.size();++row)traverse(W.view(row),d,g,W.residues[row]);}
 }
 require(R.M>=T.L,"threshold must be guaranteed nonempty or independently seeded");
 return N;
}
static void dump(const Table&T,const Positive&P,const NewResult&N,bool densecheck,double seconds){auto&R=N.result;
 std::cout<<"{\"p\":"<<T.p<<",\"r\":"<<T.r<<",\"t\":"<<T.L<<",\"B\":"<<T.B<<",\"D\":"<<T.D<<",\"maximum\":"<<R.M<<",\"prefix\":"<<R.prefix<<",\"witness\":"<<R.witness
 <<",\"seed\":"<<N.seed<<",\"seed_witness\":"<<N.seed_witness<<",\"target\":"<<N.target
 <<",\"plane_bytes\":"<<P.data.size()*4<<",\"words\":"<<N.words<<",\"stages\":"<<N.filterstats.stages<<",\"direct_plane_queries\":"<<N.filterstats.direct<<",\"signed_plane_queries\":"<<N.filterstats.signed_reads
 <<",\"prefilter_known_queries\":"<<N.filterstats.certified_queries<<",\"relaxed_bits\":"<<N.relaxed_bits<<",\"survivor_representatives\":"<<N.survivor_reps<<",\"orientations\":"<<N.orientation_total
 <<",\"head_visits\":"<<R.st.head_visits<<",\"extension_comparisons\":"<<R.st.extensions<<",\"oracle_queries\":"<<R.st.queries<<",\"euclid_steps\":"<<R.st.euclid
 <<",\"sparse_word_descriptors\":"<<N.sparse_words<<",\"nonempty_warps_for_emission\":"<<N.nonempty_warps<<",\"retained_bits_before_gcd\":"<<N.output_bits_before_gcd
 <<",\"wheel_words\":["<<N.wheel_words[1]<<","<<N.wheel_words[2]<<","<<N.wheel_words[3]<<","<<N.wheel_words[6]<<"],\"wheel_bytes\":["<<N.wheel_plane_bytes[1]<<","<<N.wheel_plane_bytes[2]<<","<<N.wheel_plane_bytes[3]<<","<<N.wheel_plane_bytes[6]<<"]"
 <<",\"dense_checked\":"<<(densecheck?"true":"false")<<",\"cpu_seconds\":"<<seconds<<",\"GPU_executed\":false}\n";
}
static void compare_dense_profile(const Table&T,const NewResult&N){
 std::vector<uint8_t> w(T.p);u64 x=1;for(u32 j=0;j<T.p-1;++j){w[x]=j%T.r;x=x*T.g%T.p;}
 u32 max=0,pre=0,witness=0;
 for(u32 a=1;a<T.p;){u32 b=a+1;while(b<T.p&&w[b]==w[a])++b;u32 l=b-a;if(l>max){max=l;witness=a;}if(a==1)pre=l;a=b;}
 require(N.result.M==max&&N.result.prefix==pre&&N.result.witness==witness,"positive dense profile");
}
#ifndef VDW_POSITIVE_NO_MAIN
int main(int argc,char**argv){try{
 if(argc==4&&(std::string(argv[1])=="--profile"||std::string(argv[1])=="--wheel")){
  u32 p=parse(argv[2]),r=parse(argv[3]);require(prime(p)&&r>=2&&r<=9&&(p-1)%r==0,"profile domain");auto sel=select_t(p,r);require(sel.t>=3,"dense fallback required");
  Table T(p,r,sel.t);Positive P(T);auto t0=std::chrono::steady_clock::now();auto N=improved(T,P,true,false,true,std::string(argv[1])=="--wheel");double s=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  if(p<=2000000)compare_dense_profile(T,N);
  u64 c=power(N.result.witness,(p-1)/r,p);for(u32 j=0;j<N.result.M;++j)require(power(N.result.witness+j,(p-1)/r,p)==c,"new witness");
  if(N.result.witness>1)require(power(N.result.witness-1,(p-1)/r,p)!=c,"left neighbour");
  if((u64)N.result.witness+N.result.M<p)require(power(N.result.witness+N.result.M,(p-1)/r,p)!=c,"right neighbour");
  dump(T,P,N,p<=2000000,s);return 0;
 }
 std::cerr<<"usage: --profile prime r\n";return 2;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
