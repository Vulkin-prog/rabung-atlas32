#define VDW_TARGET_NO_MAIN
#include "../src/target_reference.cpp"
#include <random>
struct Obs {
 const std::vector<unsigned char>* value;const std::vector<unsigned char>* available;u32 np;
 u32 valid(int bit)const{u32 out=0;for(u32 k=0;k<32;++k){int p=bit+(int)k;if(p>=0&&p<(int)value->size()&&(*available)[p])out|=1u<<k;}return out;}
 u32 window(u32 plane,int bit)const{u32 out=0;for(u32 k=0;k<32;++k){int p=bit+(int)k;if(p>=0&&p<(int)value->size()&&(((*value)[p]>>plane)&1))out|=1u<<k;}return out;}
 u32 aligned(u32 plane,u32 word)const{return window(plane,(int)(32*word));}
};
static u32 bits_for(u32 r){u32 n=0;while((1u<<n)<r)++n;return n;}
static u32 run_length(const std::vector<unsigned char>&v,u32 b,u32 step,u32 t){
 auto c=v[b];u32 l=b-(t-1)*step,h=b;
 while(l>=step&&v[l-step]==c)l-=step;
 while(h+step<v.size()&&v[h+step]==c)h+=step;
 return (h-l)/step+1;
}
static void words(){
 std::mt19937 rng(93473291);u64 masks=0,anchors=0,plane_pairs=0,partial=0;
 for(auto [c,f]:{std::pair{2u,4u},std::pair{2u,6u},std::pair{2u,8u},std::pair{3u,6u},std::pair{3u,9u},std::pair{4u,8u}}){
  u32 selected=minimal_planes(c,f),minbits=0;while((1u<<minbits)<f/c)++minbits;
  require((u32)__builtin_popcount(selected)==minbits,"attains information lower bound");
  for(u32 a=0;a<f;++a)for(u32 b=0;b<f;++b){require((a==b)==(a%c==b%c&&((encoded(a,f)^encoded(b,f))&selected)==0),"finite-fibre identity");++plane_pairs;}
  for(u32 tr=0;tr<3500;++tr){u32 n=100+rng()%100,step=1+rng()%4,t=2+rng()%6,qf=t+rng()%8,qc=qf+rng()%9;
   std::vector<unsigned char>raw(n),coarse(n),fine(n),known(n,1);for(auto&v:raw)v=rng()%f;
   u32 b=(t-1)*step+rng()%(n-(t-1)*step),colour=rng()%f;
   // Plant various coarse/fine runs; anchors can be internal and near edges.
   for(u32 j=0;j<t;++j)raw[b-j*step]=(tr&1)?colour:(colour%c+c*(rng()%(f/c)));
   for(u32 j=1;j<qc+2;++j){if(b+j*step<n&&rng()%4)raw[b+j*step]=colour;if(b>=(t-1+j)*step&&rng()%4)raw[b-(t-1+j)*step]=colour;}
   for(u32 i=0;i<n;++i){coarse[i]=encoded(raw[i]%c,c);fine[i]=encoded(raw[i],f);if(tr%3==0)known[i]=rng()%3!=0;}
   for(u32 j=0;j<t;++j)known[b-j*step]=1;
   Obs C{&coarse,&known,bits_for(c)},F{&fine,&known,bits_for(f)};
   for(u32 w=0;32*w<n;++w){u32 mask=0;
    for(u32 k=0;k<32&&32*w+k<n;++k){u32 x=32*w+k;if(x<(t-1)*step)continue;bool ok=true;for(u32 j=0;j<t;++j)if(coarse[x-j*step]!=coarse[x]||!known[x-j*step])ok=false;if(ok)mask|=1u<<k;}
    auto nest=vdw_target::two_thresholds(C,step,t,qf,qc,w,mask);
    require(nest.fine==vdw_cut::screen(C,step,t,qf,w,mask,1),"weak mask equals separate pass");
    require(nest.coarse==vdw_cut::screen(C,step,t,qc,w,mask,1),"strong mask equals separate pass");
    auto out=vdw_target::split(C,F,selected,step,t,qc,qf,w,mask);
    require(!(out.coarse&out.fine)&&!((out.coarse|out.fine)&~mask),"partition");
    for(u32 k=0;k<32;++k)if(mask&(1u<<k)){
     u32 x=32*w+k;bool same=true;for(u32 j=0;j<t;++j)if(fine[x-j*step]!=fine[x])same=false;
     if(run_length(coarse,x,step,t)>=qc)require(out.coarse&(1u<<k),"coarse witness retained");
     if(same&&run_length(fine,x,step,t)>=qf)require((out.coarse|out.fine)&(1u<<k),"fine witness retained");
     if(out.fine&(1u<<k))require(same,"fine-only path has fine anchor");
     ++anchors;if(tr%3==0)++partial;
    }++masks;
   }
  }
 }
 // A fine-only survivor need not have a coarse-surviving mask: removing that
 // branch loses valid fine maxima. Entire known window has length 5<qc=8.
 std::vector<unsigned char>raw={3,1,0,0,0,0,0,1,3},c(raw.size()),f(raw.size()),known(raw.size(),1);
 for(u32 i=0;i<raw.size();++i){c[i]=raw[i]%2;f[i]=raw[i];}
 Obs C{&c,&known,1},F{&f,&known,2};auto out=vdw_target::split(C,F,2,1,3,8,5,0,1u<<4);
 require(out.coarse==0&&out.fine==(1u<<4),"negative coarse-only failure");
 // The additional bit is necessary: labels 0 and 2 agree modulo 2.
 require(((encoded(0,4)^encoded(2,4))&minimal_planes(2,4))!=0,"missing residual bit rejects falsely fine anchor");
 std::cout<<"{\"status\":\"PASS\",\"row_masks\":"<<masks<<",\"anchor_checks\":"<<anchors<<",\"partial_observation_anchors\":"<<partial<<",\"finite_colour_pairs\":"<<plane_pairs<<",\"negative_cases\":2}\n";
}
static auto signature(const SharedOutput&O){return std::make_tuple(O.coarse_result.maximum,O.coarse_result.prefix,O.coarse_result.witness,O.fine_result.maximum,O.fine_result.prefix,O.fine_result.witness);}
static void profiles(std::string suite){std::vector<u32>ps;
 if(suite=="small"){for(u32 p=101;p<1000;++p)if(prime(p))ps.push_back(p);}
 else if(suite=="million"){for(u32 p=1000000;p<1002000&&ps.size()<12;++p)if(prime(p))ps.push_back(p);}
 else throw std::invalid_argument("suite");
 u64 n=0,fallbacks=0;
 for(auto p:ps)for(auto[c,f]:{std::pair{2u,4u},std::pair{2u,6u},std::pair{2u,8u},std::pair{3u,6u},std::pair{3u,9u},std::pair{4u,8u}})if((p-1)%f==0){
  try{Prepared F(p,f,16);(void)F;}catch(const std::exception&e){if(std::string(e.what())=="fallback required"){++fallbacks;continue;}throw;}
  vdw_cut::Counts cut;auto base=shared_screened(p,c,f,1,16,3,64,false,1,&cut);
  for(u32 t:{0u,1u})for(u32 batch:{1u,64u,1024u}){TargetCounts stats;auto O=targeted_run(p,c,f,t,16,3,batch,stats);require(signature(base)==signature(O),"baseline exact equality");++n;}
 }
 std::cout<<"{\"status\":\"PASS\",\"suite\":\""<<suite<<"\",\"coupled_runs\":"<<n<<",\"individual_profiles\":"<<2*n<<",\"explicit_fallbacks\":"<<fallbacks<<",\"dense_compared\":true}\n";
}
int main(int argc,char**argv){try{std::string mode=argc>1?argv[1]:"words";if(mode=="words")words();else profiles(mode);}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
