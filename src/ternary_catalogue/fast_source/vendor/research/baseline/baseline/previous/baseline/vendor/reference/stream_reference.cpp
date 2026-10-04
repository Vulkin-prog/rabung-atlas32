// CPU research reference. The retained prior Table is also used as a full-plane
// oracle; the proposed GPU filter stores ONLY the positive folded planes.
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmisleading-indentation"
#endif
#define main retained_prior_main
#include "prior_word32.cpp"
#undef main
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#include "folded_filter.hpp"
#include <map>
#include <random>
#include <set>
#include <sstream>
struct Fold {
 const Table& T; u32 stride,np,two,v2,flip; std::vector<u32> data;
 explicit Fold(const Table& t):T(t),stride((t.B+32)/32+2),np(t.planes.size()),two(1),v2(0),flip(0) {
  while((T.r%(two*2))==0){two*=2;++v2;}
  if(T.sigma) flip=two/2;
  data.assign(np*stride,0);
  for(u32 x=1;x<=T.B;++x){u32 code=encode(T.c[x]);for(u32 j=0;j<np;++j)if(code&(1u<<j))data[j*stride+(x>>5)]|=1u<<(x&31);}
 }
 u32 encode(u32 c)const{return (c%two)+two*(c%(T.r/two));}
 vdw::FoldView view()const{return {data.data(),stride,T.B,np,flip};}
};
struct Selection {u32 t=1;u64 lower=0,upper=0;};
static Selection select_t(u32 p,u32 r){u64 H=(u64)std::sqrt((long double)p);while(H*H<p)++H;
 Selection S;u64 Q=1;for(u32 t=2;t<256&&t<p;++t){Q*=r;__uint128_t E=(__uint128_t)(Q-1)*(t-1)*H;if(E>=p-t)break;
 S={t,(p-t-(u64)E+Q-1)/Q,(p-t+(u64)E)/Q};}return S;}
struct Stats {u64 words=0,stages=0,raw=0,primitive=0,oriented=0,head_visits=0,queries=0,euclid=0,short_hits=0,extensions=0;};
// Arbitrary nonzero residue -> colour using the SMALL table, not powmod.
static int small_oracle(const Table&T,u32 x,Stats&st){require(x>0&&x<T.p,"oracle domain");++st.queries;
 if(x<=T.B){++st.short_hits;return T.c[x];}
 if(T.p-x<=T.B){++st.short_hits;return (T.c[T.p-x]+T.sigma)%T.r;}
 u64 a=T.p,b=x;i64 c=0,d=1;
 while(b>T.B){u64 q=a/b,rem=a-q*b;i64 next=c-(i64)q*d;a=b;b=rem;c=d;d=next;++st.euclid;}
 u64 ad=d<0?(u64)-d:(u64)d;require(b&&ad&&b<=T.B&&ad<=T.B,"small quotient bound");
 int bc=(T.c[ad]+(d<0?T.sigma:0))%T.r;return (T.c[b]+T.r-bc)%T.r;
}
struct StreamResult{u32 M=0,prefix=0,witness=0;Stats st;u64 multiplicity_bound=0,upper=0,audit_max_mult=0;std::map<u32,u32> audit_heads;std::vector<u32> audit_starts;};
static u64 mult_bound(const Table&T){u64 h=T.L-1,q=(4ull*T.B*T.B)/(h*T.p)+1;return (q*q-1)/2;}
static StreamResult streaming(const Table&T,const Fold&F,bool audit=false,bool compare=false){
 StreamResult R;R.multiplicity_bound=mult_bound(T);R.upper=select_t(T.p,T.r).upper;
 std::vector<u32> inv(T.D+1);inv[1]=1;for(u32 d=2;d<=T.D;++d)inv[d]=T.p-(u64)(T.p/d)*inv[T.p%d]%T.p;
 auto visit=[&](u32 x,i64 a,u32 d){require(x>=1&&(u64)x+T.L<=T.p,"linear lift");++R.st.oriented;
  if(audit)R.audit_starts.push_back(x);
  const int col=(T.colour(a)+T.r-T.c[d])%T.r;
  // Only maximal-run HEADS are completed. Duplicates are harmless for max.
  if(x>1){i64 pre=a-d;int cp;
   if(pre>=-(i64)T.B&&pre<=T.B&&pre){cp=(T.colour(pre)+T.r-T.c[d])%T.r;}
   else cp=small_oracle(T,x-1,R.st);
   if(cp==col)return;
  }
  ++R.st.head_visits;u32 len=T.L;
  while((u64)x+len<T.p){i64 z=a+(i64)len*d;int cc;
   if(z&&z>=-(i64)T.B&&z<=T.B)cc=(T.colour(z)+T.r-T.c[d])%T.r;
   else cc=small_oracle(T,x+len,R.st);
   ++R.st.extensions;if(cc!=col)break;++len;
  }
  if(len>R.M||(len==R.M&&x<R.witness)){R.M=len;R.witness=x;}
  if(x==1)R.prefix=len;
  if(audit){auto it=R.audit_heads.find(x);require(it==R.audit_heads.end()||it->second==len,"head length consistent");R.audit_heads[x]=len;}
 };
 const auto order=shifts(T.L,false);
 for(u32 d=1;d<=T.D;++d){i64 last=2*(i64)T.B-(T.L-1)*(i64)d;if(last<0)continue;
  u32 nw=(last/2+32)/32;for(u32 w=0;w<nw;++w){++R.st.words;unsigned long long stages=0;u32 bits=vdw::folded_candidates(F.view(),d,T.L,w,&stages);R.st.stages+=stages;
   if(compare){u64 dummy=0;auto expected=candidates(T,d,T.L,w,order,dummy); if(bits!=expected){std::cerr<<"p="<<T.p<<" r="<<T.r<<" t="<<T.L<<" B="<<T.B<<" d="<<d<<" w="<<w<<" got="<<bits<<" want="<<expected<<"\n";} require(bits==expected,"folded/full masks");}
   while(bits){u32 k=__builtin_ctz(bits);bits&=bits-1;++R.st.raw;i64 a=32*(i64)w+k-T.B;
    if(std::gcd((u64)(a<0?-a:a),(u64)d)!=1)continue;
    ++R.st.primitive;u64 aa=(a%(i64)T.p+T.p)%T.p;u32 x=aa*inv[d]%T.p;
    visit(x,a,d);u32 xr=T.p-x-(T.L-1);if(xr!=x)visit(xr,-a-(T.L-1)*(i64)d,d);
   }
  }
 }
 require(R.M>=T.L,"no long run: caller must use dense fallback / valid threshold");
 if(!R.prefix){while(R.prefix<T.B&&T.c[R.prefix+1]==0)++R.prefix;require(R.prefix<T.L,"missing prefix head");}
 if(audit){std::sort(R.audit_starts.begin(),R.audit_starts.end());u64 maxi=0;
  for(size_t a=0;a<R.audit_starts.size();){size_t b=a+1;while(b<R.audit_starts.size()&&R.audit_starts[b]==R.audit_starts[a])++b;maxi=std::max<u64>(maxi,b-a);a=b;}
  require(maxi<=R.multiplicity_bound,"primitive multiplicity bound");R.audit_max_mult=maxi;
  R.audit_starts.erase(std::unique(R.audit_starts.begin(),R.audit_starts.end()),R.audit_starts.end());
 }
 require(R.st.oriented<=R.multiplicity_bound*R.upper,"visit upper bound");
 require(R.st.extensions<=R.multiplicity_bound*R.upper,"extension upper bound");
 return R;
}
static u64 check_dense(const Table&T,const Fold&F,const StreamResult&R){std::vector<uint8_t> dense(T.p);u64 x=1;
 for(u32 j=0;j<T.p-1;++j){dense[x]=j%T.r;x=x*T.g%T.p;}
 for(u32 a=1;a<=T.B;++a)require(T.c[a]==dense[a],"table against multiplicative walk");
 u64 window_checks=0;auto view=F.view();
 for(i64 a=-(i64)T.B;a<=(i64)T.B;a+=1)for(u32 plane=0;plane<F.np;++plane){u32 word=view.window(plane,(int)a);
  for(u32 j=0;j<32;++j){i64 z=a+j;if(!z||z<-(i64)T.B||z>T.B)continue;
   u32 want=(F.encode(T.colour(z))>>plane)&1;require(((word>>j)&1u)==want,"folded bit reconstruction");++window_checks;}}
 std::vector<u32> starts;std::map<u32,u32> heads;u32 M=0,prefix=0;
 for(u32 a=1;a<T.p;){u32 b=a+1;while(b<T.p&&dense[b]==dense[a])++b;u32 len=b-a;M=std::max(M,len);if(a==1)prefix=len;
  if(len>=T.L){heads[a]=len;for(u32 q=a;(u64)q+T.L<=b;++q)starts.push_back(q);}a=b;}
 require(R.M==M&&R.prefix==prefix&&R.audit_heads==heads&&R.audit_starts==starts,"streaming/full lists or heads");return window_checks;
}
static void print_profile(const Table&T,const Fold&F,const StreamResult&R,bool dense_checked){
 std::cout<<"{\"p\":"<<T.p<<",\"r\":"<<T.r<<",\"t\":"<<T.L<<",\"B\":"<<T.B<<",\"D\":"<<T.D<<",\"maximum\":"<<R.M<<",\"prefix\":"<<R.prefix<<",\"witness\":"<<R.witness
 <<",\"folded_plane_bytes\":"<<F.data.size()*4<<",\"signed_plane_bytes\":"<<T.planes.size()*T.planes[0].size()*4
 <<",\"words32\":"<<R.st.words<<",\"stages32\":"<<R.st.stages<<",\"raw_hits\":"<<R.st.raw<<",\"primitive_hits\":"<<R.st.primitive<<",\"oriented_visits\":"<<R.st.oriented<<",\"head_visits_with_duplicates\":"<<R.st.head_visits
 <<",\"extension_comparisons\":"<<R.st.extensions<<",\"oracle_queries\":"<<R.st.queries<<",\"oracle_euclid_steps\":"<<R.st.euclid<<",\"oracle_direct_small_hits\":"<<R.st.short_hits
 <<",\"audit_unique_starts\":"<<R.audit_starts.size()<<",\"audit_unique_heads\":"<<R.audit_heads.size()<<",\"audit_max_multiplicity\":"<<R.audit_max_mult<<",\"primitive_multiplicity_upper\":"<<R.multiplicity_bound<<",\"dense_oracle\":"<<(dense_checked?"true":"false")<<",\"production_start_sort\":false,\"GPU_executed\":false}"<<'\n';
}
static u32 parse(const char*s){std::string a=s;require(!a.empty()&&a.find_first_not_of("0123456789")==std::string::npos,"unsigned decimal required");size_t n=0;u64 x=std::stoull(a,&n);require(n==a.size()&&x<=UINT32_MAX,"u32 overflow");return x;}
int stream_entry(int argc,char**argv){try{
 if(argc==2&&std::string(argv[1])=="--test"){
  u64 profiles=0,windows=0,oracle_checks=0;
  for(u32 p=5;p<=1000;++p)if(prime(p))for(u32 r=2;r<=9;++r)if((p-1)%r==0){auto sel=select_t(p,r);if(sel.t<3)continue;Table T(p,r,sel.t);Fold F(T);auto R=streaming(T,F,true,true);windows+=check_dense(T,F,R);Stats st;
   for(u32 x=1;x<p;++x){u64 y=power(x,(p-1)/r,p);int c=small_oracle(T,x,st);require(y==power(T.g,(p-1)/r*c,p),"independent pow colour");++oracle_checks;}++profiles;}
  std::cout<<"{\"status\":\"PASS\",\"profiles\":"<<profiles<<",\"folded_bit_comparisons\":"<<windows<<",\"oracle_pow_comparisons\":"<<oracle_checks<<",\"GPU_executed\":false}\n";return 0;
 }
 if((argc==4||argc==6)&&std::string(argv[1])=="--profile"){
  u32 p=parse(argv[2]),r=parse(argv[3]);require(prime(p)&&r>=2&&r<=9&&(p-1)%r==0,"profile domain");auto s=select_t(p,r);require(s.t>=3,"use dense fallback for small input");Table T(p,r,s.t);Fold F(T);bool audit=p<=2000000||argc==6;require(argc!=6||std::string(argv[4])=="--audit-list","expected --audit-list path");auto R=streaming(T,F,audit,p<=2000000);if(p<=2000000)check_dense(T,F,R);
  if(argc==6){std::ofstream out(argv[5],std::ios::binary);require((bool)out,"audit output open");for(u32 x:R.audit_starts){char b[4];for(u32 j=0;j<4;++j)b[j]=(x>>(8*j))&255;out.write(b,4);}require((bool)out,"audit output write");}
  u64 z=power(R.witness,(p-1)/r,p);for(u32 j=0;j<R.M;++j)require(power(R.witness+j,(p-1)/r,p)==z,"maximum witness pow");
  if(R.witness>1)require(power(R.witness-1,(p-1)/r,p)!=z,"left witness");
  if((u64)R.witness+R.M<p)require(power(R.witness+R.M,(p-1)/r,p)!=z,"right witness");
  print_profile(T,F,R,p<=2000000);return 0;
 }
 std::cerr<<"Usage: --test | --profile prime r\n";return 2;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}

#ifndef VDW_STREAM_NO_MAIN
int main(int argc,char**argv){return stream_entry(argc,argv);}
#endif
