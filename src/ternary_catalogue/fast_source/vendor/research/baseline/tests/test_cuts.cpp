#define VDW_CUT_NO_MAIN
#include "../src/cut_reference.cpp"
#include <random>
// Independent observation model. No rational table and no folded reader here.
struct Obs {
 const std::vector<unsigned char>* value;const std::vector<unsigned char>* available;u32 np;
 u32 valid(int bit)const{u32 out=0;for(u32 k=0;k<32;++k){int pos=bit+(int)k;
  if(pos>=0&&pos<(int)value->size()&&(*available)[pos])out|=1u<<k;}return out;}
 u32 window(u32 plane,int bit)const{u32 out=0;for(u32 k=0;k<32;++k){int pos=bit+(int)k;
  if(pos>=0&&pos<(int)value->size()&&(((*value)[pos]>>plane)&1))out|=1u<<k;}return out;}
 u32 aligned(u32 plane,u32 word)const{return window(plane,32*word);}
};
static bool scalar_possible(const Obs&v,u32 b,u32 step,u32 t,u32 q){
 unsigned colour=(*v.value)[b];
 // Examine every placement of a q-window around the t-anchor. Unknowns,
 // including positions outside this finite test word, are not mismatches.
 for(u32 left=0;left<=q-t;++left){bool possible=true;
  for(int j=-(int)left;j<(int)(q-left);++j){int pos=(int)b-j*(int)step;
   if(pos>=0&&pos<(int)v.value->size()&&(*v.available)[pos]&&(*v.value)[pos]!=colour){possible=false;break;}}
  if(possible)return true;
 }
 return false;
}
static u64 mask_tests(){std::mt19937 rng(3907729);u64 checked=0;
 for(u32 alph: {2u,3u,4u,6u,8u,9u}){
  u32 np=0;while((1u<<np)<alph)++np;
  for(u32 trial=0;trial<300;++trial){u32 n=65+rng()%150,step=1+rng()%5,t=2+rng()%7,q=t+1+rng()%10;
   std::vector<unsigned char>w(n),known(n,1);for(auto&c:w)c=rng()%alph;
   // Plant an anchor, sometimes embedded in a much longer run.
   u32 b=(t-1)*step+rng()%(n-(t-1)*step);unsigned c=rng()%alph;
   for(u32 j=0;j<t;++j)w[b-j*step]=c;
   for(u32 j=1;j<=q;++j){if(b+j*step<n&&rng()%3)w[b+j*step]=c;if(b>=(t-1+j)*step&&rng()%3)w[b-(t-1+j)*step]=c;}
   for(u32 i=0;i<n;++i)known[i]=rng()%4!=0;for(u32 j=0;j<t;++j)known[b-j*step]=1;
   Obs v{&w,&known,np};
   for(u32 word=0;word*32<n;++word){u32 mask=0;
    for(u32 k=0;k<32&&32*word+k<n;++k){u32 x=32*word+k;if(x<(t-1)*step)continue;bool ok=true;
     for(u32 j=0;j<t;++j)if(w[x-j*step]!=w[x]||!known[x-j*step])ok=false;if(ok)mask|=1u<<k;}
    auto exact=vdw_cut::screen(v,step,t,q,word,mask,4);u32 want=0;
    for(u32 k=0;k<32;++k)if(mask&(1u<<k)){if(scalar_possible(v,word*32+k,step,t,q))want|=1u<<k;++checked;}
    require(exact==want,"exact observation mask/scalar windows");
    for(u32 mode: {1u,2u,3u}){auto got=vdw_cut::screen(v,step,t,q,word,mask,mode);require((want&~got)==0,"no false negative");require((got&~mask)==0,"no invented anchor");}
   }
  }
 }
 return checked;
}
static u64 lemma_tests(){u64 checked=0;
 // Enumerate all binary words up to 13; examine every constant anchor t>=2
 // and every attainable surrounding pair of cut distances.
 for(u32 n=2;n<=13;++n)for(u32 code=0;code<(1u<<n);++code){
  std::vector<unsigned char>w(n);for(u32 i=0;i<n;++i)w[i]=(code>>i)&1;
  for(u32 a=0;a<n;++a){u32 h=a;while(h<n&&w[h]==w[a])++h;u32 l=a;while(l&&w[l-1]==w[a])--l;
   for(u32 t=2;t<=h-a;++t)for(u32 q=t+1;q<=n;++q){u32 U=(q-t+1)/2,V=q-t+1-U;
    bool left=false,right=false;
    for(u32 i=1;i<=U&&i<=a;++i)left|=w[a-i]!=w[a];
    for(u32 i=1;i<=V&&a+t-1+i<n;++i)right|=w[a+t-1+i]!=w[a];
    if(left&&right)require(h-l<q,"bilateral cuts lemma");++checked;
   }
  }
 }
 return checked;
}
static void negative_tests(){
 // Either ONE cut alone or replacing unknowns with differences would reject
 // a genuine long run: these are guards on the intended logic.
 std::vector<unsigned char>w={1,0,0,0,0,0,1},known(w.size(),1);
 Obs v{&w,&known,1};u32 anchor=1u<<3; // t=3: positions 3,2,1, real run length 5
 require(vdw_cut::screen(v,1,3,5,0,anchor,1)==anchor,"one cut is insufficient");
 known[0]=known[4]=known[5]=known[6]=0;
 require(vdw_cut::screen(v,1,3,5,0,anchor,4)==anchor,"unknown positions cannot reject");
 require(vdw_cut::screen(v,1,3,3,0,anchor,1)==anchor,"equal threshold needs no pruning");
 // A length q-1 run may be rejected, whereas q must survive (strictness).
 w={1,1,0,0,0,0,1};known.assign(w.size(),1);Obs z{&w,&known,1};anchor=1u<<5;
 require(vdw_cut::screen(z,1,3,5,0,anchor,1)==0,"short run should reject");
 require(vdw_cut::screen(z,1,3,4,0,anchor,1)==anchor,"tie length must survive");
}
static void profiles(std::string mode){require(mode=="small"||mode=="million","unknown suite");u64 done=0,fall=0;std::vector<u32>ps;
 if(mode=="small"){for(u32 p=101;p<1300;++p)if(prime(p))ps.push_back(p);}
 else{for(u32 p=1000000;p<1001000&&ps.size()<12;++p)if(prime(p))ps.push_back(p);}
 for(auto p:ps)for(auto [c,f]:{std::pair{2u,4u},std::pair{2u,6u},std::pair{2u,8u},std::pair{4u,8u},std::pair{3u,6u},std::pair{3u,9u}})if((p-1)%f==0){
  try{Prepared probe(p,f,16);(void)probe;}catch(const std::exception&e){if(std::string(e.what())=="fallback required"){++fall;continue;}throw;}
  auto want=shared_run(p,c,f,0,16,3,64,true);
  for(u32 tm: {0u,1u})for(u32 scr: {1u,4u})for(u32 batch: {1u,64u}){
   vdw_cut::Counts st;auto got=shared_screened(p,c,f,tm,16,3,batch,true,scr,&st);
   auto sig=[](const SharedOutput&x){return std::make_tuple(x.coarse_result.maximum,x.coarse_result.prefix,x.coarse_result.witness,x.fine_result.maximum,x.fine_result.prefix,x.fine_result.witness);};
   require(sig(got)==sig(want),"screened/baseline profile mismatch");++done;
  }
 }
 std::cout<<"{\"status\":\"PASS\",\"suite\":\""<<mode<<"\",\"coupled_profiles\":"<<done<<",\"individual_profiles\":"<<2*done<<",\"seed_fallbacks\":"<<fall<<",\"dense_oracle\":true,\"cuda_compiled\":false,\"gpu_executed\":false}\n";
}
int main(int argc,char**argv){try{std::string mode=argc>1?argv[1]:"words";
 if(mode=="words"){auto a=lemma_tests(),b=mask_tests();negative_tests();std::cout<<"{\"status\":\"PASS\",\"binary_word_anchor_threshold_cases\":"<<a<<",\"partial_observation_anchors\":"<<b<<",\"negative_cases\":5}\n";}
 else profiles(mode);
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
