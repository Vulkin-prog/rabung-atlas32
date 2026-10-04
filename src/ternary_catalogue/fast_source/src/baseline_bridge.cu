// The source below is extracted READ-ONLY from the user's local repository and
// accepted only by its exact Git blob ID. No edits to that source are performed.
#define main archived_scanner_main
#include "../build/original_scan_gpu.cu"
#undef main
static Buffers kit_buffers{};static bool kit_ready=false;
extern "C" void baseline_init(unsigned cap){if(!kit_ready){alloc_buffers_b7(kit_buffers,(long)(cap/2ull+128));kit_ready=true;}}
template<unsigned MASK>static void subset(unsigned p,unsigned g,int*ans){
 auto &b=kit_buffers;u32 nh=(p-1)/2,m=1;
 for(int r=2;r<=9;++r)if(MASK&(1u<<(r-2)))m=h_lcm(m,r);
 u32 e=(p-1)/m,ni=mont_ninv(p),r2=mont_r2(p),one=(u32)((1ull<<32)%p);
 u32 z=h_powmod(g,e,p),cur=1;std::vector<std::pair<u32,u32>> tab(m);
 for(u32 j=0;j<m;++j){tab[j]={montmul(cur,r2,p,ni),j};cur=h_mulmod(cur,z,p);}std::sort(tab.begin(),tab.end());
 std::vector<u32> val(m),lab(m);for(u32 j=0;j<m;++j){val[j]=tab[j].first;lab[j]=tab[j].second;}
 CK(cudaMemcpy(b.d_sval,val.data(),m*4,cudaMemcpyHostToDevice));CK(cudaMemcpy(b.d_sj,lab.data(),m*4,cudaMemcpyHostToDevice));
 u32 nseg=(nh+SEG-1)/SEG;
 k_scan_char_best<MASK><<<(nseg+127)/128,128,2*(size_t)m*4>>>(p,ni,r2,one,e,m,b.d_sval,b.d_sj,nh,nseg,b.summ);
 ordered_reduce_active(b,nseg,MASK);CK(cudaMemcpy(b.h_sout,b.d_sout,NR*sizeof(Summ),cudaMemcpyDeviceToHost));
 for(int r=2;r<=9;++r)if(MASK&(1u<<(r-2))){auto s=b.h_sout[r-2];ans[r]=(nh%r==0)?std::max(s.best,2*s.sL):s.best;}
}
extern "C" int baseline_profile(unsigned p,unsigned mask,unsigned*out,unsigned*prefix,int full){
 int a[10]{},m[10]{};u32 g=primitive_root(p);
 if(full)process_prime_char_v2a(kit_buffers,p,g,m,a,nullptr);
 else switch(mask){case 1:subset<1>(p,g,m);break;case 5:subset<5>(p,g,m);break;case 18:subset<18>(p,g,m);break;case 19:subset<19>(p,g,m);break;case 23:subset<23>(p,g,m);break;default:return 2;}
 for(unsigned r=2;r<=9;++r)if(mask&(1u<<(r-2))){out[r]=m[r];u32 k=1;while(k+1<p&&h_powmod(k+1,(p-1)/r,p)==1)++k;prefix[r]=k;}
 return 0;
}
