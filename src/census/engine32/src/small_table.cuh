#pragma once
#include <cuda_runtime.h>
#include "completion_core.hpp"
#include "mont32.hpp"
namespace kit {
__device__ __forceinline__ u32 modpow(u32 x,u32 e,u32 p,u32 ni,u32 one){if(!e)return one;u32 z=x;int b=31-__clz(e);for(--b;b>=0;--b){z=mont(z,z,p,ni);if((e>>b)&1)z=mont(z,x,p,ni);}return z;}
__global__ void small_table(u32 p,u32 B,u32 m,u32 ni,u32 R2,u32 one,const u32*roots,const u32*labels,unsigned short*out,unsigned*failed){
 u32 x=blockIdx.x*blockDim.x+threadIdx.x;if(x>B)return;if(!x){out[0]=0;return;}
 u32 y=modpow(mont(x,R2,p,ni),(p-1)/m,p,ni,one),lo=0,hi=m;
 while(lo<hi){u32 mid=(lo+hi)/2;if(roots[mid]<y)lo=mid+1;else hi=mid;}
 if(lo==m||roots[lo]!=y){atomicExch(failed,1u);out[x]=0;}else out[x]=(unsigned short)labels[lo];
}
__global__ void inverse_table(u32 p,u32 D,u32 ni,u32 R2,u32 one,u32*out){u32 d=blockIdx.x*blockDim.x+threadIdx.x;if(d>D)return;if(!d){out[d]=0;return;}u32 v=modpow(mont(d,R2,p,ni),p-2,p,ni,one);
#if VDW_COMPLETION_MONT
 out[d]=v;
#else
 out[d]=mont(v,1,p,ni);
#endif
}
__global__ void fine_labels(u32 B,u32 f,const unsigned short*base,unsigned char*out){u32 x=blockIdx.x*blockDim.x+threadIdx.x;if(x<=B)out[x]=base[x]%f;}
__global__ void pack_planes(u32 B,u32 r,u32 g,u32 np,u32 stride,const unsigned char*colour,u32*out){
 u32 k=blockIdx.x*blockDim.x+threadIdx.x,nr=(g==1||g==2)?1:2,total=nr*np*stride;if(k>=total)return;
 u32 word=k%stride,plane=(k/stride)%np,row=k/(stride*np),res=g==1?0:row==0?1:g-1,mask=0,two=1;
 while(r%(two*2)==0)two*=2;
 if(word)for(u32 bit=0;bit<32;++bit){u64 x=(u64)g*((word-1)*32+bit)+res;if(x&&x<=B){u32 c=colour[x]%r,e=c%two+two*(c%(r/two));if(e&(1u<<plane))mask|=1u<<bit;}}
 out[k]=mask;
}
#ifndef KIT_WARP_PACK
#define KIT_WARP_PACK 0
#endif
// One warp reads a contiguous 32-entry word, then ballots every coarse/fine
// plane. Current groups have fine=coarse or fine=2*coarse.
__global__ void pack_pair(u32 B,u32 coarse,u32 fine,u32 g,u32 np,u32 nf,u32 stride,
 const unsigned char*colour,u32*cp,u32*fp){
 const u32 lane=threadIdx.x&31,j=(blockIdx.x*blockDim.x+threadIdx.x)/32;
 const u32 rows=(g==1||g==2)?1:2;
 if(j>=rows*stride)return;
 u32 word=j%stride,row=j/stride,res=g==1?0:row==0?1:g-1;
 u64 x=word?(u64)g*((word-1)*32ull+lane)+res:0;
 bool valid=word&&x&&x<=B;
 u32 raw=valid?colour[x]:0,co=raw>=coarse?raw-coarse:raw;
 u32 ec=coarse==6?(co&1)+2*(co>=3?co-3:co):co;
 u32 ef=fine==6?(raw&1)+2*(raw>=3?raw-3:raw):raw;
 for(u32 bit=0;bit<np;++bit){u32 mask=__ballot_sync(FULL,valid&&((ec>>bit)&1));if(!lane)cp[row*np*stride+bit*stride+word]=mask;}
#if KIT_SINGLETON_FAST
 if(coarse==fine)return; // fp aliases cp for singleton groups.
#endif
#if KIT_COMPACT_FINE
 if(coarse!=fine){u32 residual=coarse==2?1:0;u32 mask=__ballot_sync(FULL,valid&&((ef>>residual)&1));if(!lane)fp[row*stride+word]=mask;return;}
#endif
 for(u32 bit=0;bit<nf;++bit){u32 mask=__ballot_sync(FULL,valid&&((ef>>bit)&1));if(!lane)fp[row*nf*stride+bit*stride+word]=mask;}
}
#ifndef KIT_WARP_SEED
#define KIT_WARP_SEED 0
#endif
#ifndef KIT_SINGLETON_SEED
#define KIT_SINGLETON_SEED 0
#endif
#if KIT_SINGLETON_SEED && !KIT_WARP_SEED
#error Singleton seed requires warp seed
#endif
#if KIT_SINGLETON_SEED
// Equal coarse/fine profiles share the same segment maximum and prefix.
__global__ void seed_single(u32 p,u32 B,const unsigned char*col,State*s){
 const u32 lane=threadIdx.x&31,warp=(blockIdx.x*blockDim.x+threadIdx.x)/32;
 u32 a=1+warp*1024;if(a>B)return;
 u32 carry=0,last=~0u;u64 best=0;
 for(u32 tile=0;tile<32;++tile){
  u32 x=a+tile*32+lane;bool valid=x<=B;
  u32 v=valid?col[x]:~0u,prev=__shfl_up_sync(FULL,v,1);if(!lane)prev=last;
  u32 changes=__ballot_sync(FULL,v!=prev),mask=changes&(FULL>>(31-lane));
  u32 run=mask?lane-(31-__clz(mask))+1:carry+lane+1;
  if(valid){u32 st=x-run+1;u64 mirror=(u64)p-x;if(mirror<st)st=(u32)mirror;
   u64 z=vdw_next::pack_best(run,st);if(z>best)best=z;
  }
  carry=__shfl_sync(FULL,run,31);last=__shfl_sync(FULL,v,31);
 }
 u32 high=__reduce_max_sync(FULL,(u32)(best>>32));
 u32 low=__reduce_max_sync(FULL,(u32)(best>>32)==high?(u32)best:0u);
 if(!lane)atomicMax(&s->best_f,((unsigned long long)high<<32)|low);
}
#endif
#if KIT_WARP_SEED
__global__ void seed_runs(u32 p,u32 B,u32 c,u32 f,const unsigned char*col,State*s){
 const u32 lane=threadIdx.x&31,warp=(blockIdx.x*blockDim.x+threadIdx.x)/32,nw=gridDim.x*(blockDim.x/32);
 u32 j=warp,a=1+j*1024;if(a>B)return;
 u32 carry_c=0,carry_f=0,last_c=~0u,last_f=~0u;u64 bc=0,bf=0;
 for(u32 tile=0;tile<32;++tile){
  u32 x=a+tile*32+lane;bool valid=x<=B;
  u32 v=valid?col[x]:~0u,co=valid?(v>=c?v-c:v):~0u;
  u32 prev_c=__shfl_up_sync(FULL,co,1),prev_f=__shfl_up_sync(FULL,v,1);
  if(!lane){prev_c=last_c;prev_f=last_f;}
  u32 changes_c=__ballot_sync(FULL,co!=prev_c),changes_f=__ballot_sync(FULL,v!=prev_f);
  u32 mask_c=changes_c&(FULL>>(31-lane)),mask_f=changes_f&(FULL>>(31-lane));
  u32 rc=mask_c?lane-(31-__clz(mask_c))+1:carry_c+lane+1;
  u32 rf=mask_f?lane-(31-__clz(mask_f))+1:carry_f+lane+1;
  if(valid){
   u32 stc=x-rc+1,stf=x-rf+1,mirror=p-x;
   if(mirror<stc)stc=mirror;if(mirror<stf)stf=mirror;
   u64 z=vdw_next::pack_best(rc,stc);if(z>bc)bc=z;
   z=vdw_next::pack_best(rf,stf);if(z>bf)bf=z;
  }
  carry_c=__shfl_sync(FULL,rc,31);carry_f=__shfl_sync(FULL,rf,31);
  last_c=__shfl_sync(FULL,co,31);last_f=__shfl_sync(FULL,v,31);
 }
 u32 hc=__reduce_max_sync(FULL,(u32)(bc>>32));
 u32 lc=__reduce_max_sync(FULL,(u32)(bc>>32)==hc?(u32)bc:0u);
 u32 hf=__reduce_max_sync(FULL,(u32)(bf>>32));
 u32 lf=__reduce_max_sync(FULL,(u32)(bf>>32)==hf?(u32)bf:0u);
 if(!lane){
  atomicMax(&s->best_c,((unsigned long long)hc<<32)|lc);
  atomicMax(&s->best_f,((unsigned long long)hf<<32)|lf);
 }
}
#else
__global__ void seed_runs(u32 p,u32 B,u32 c,u32 f,const unsigned char*col,State*s){
 u32 j=blockIdx.x*blockDim.x+threadIdx.x,a=1+j*1024;if(a>B)return;u32 hi=a+1024;if(hi>B+1)hi=B+1;
 u32 rc=0,rf=0,prev_c=~0u,prev_f=~0u;u64 bc=0,bf=0;
 for(u32 x=a;x<hi;++x){u32 v=col[x],co=v%c;rc=co==prev_c?rc+1:1;rf=v==prev_f?rf+1:1;prev_c=co;prev_f=v;
  u32 stc=x-rc+1,stf=x-rf+1,mirror=p-x;if(mirror<stc)stc=mirror;if(mirror<stf)stf=mirror;
  u64 z=vdw_next::pack_best(rc,stc);if(z>bc)bc=z;z=vdw_next::pack_best(rf,stf);if(z>bf)bf=z;}
 atomicMax(&s->best_c,(unsigned long long)bc);atomicMax(&s->best_f,(unsigned long long)bf);
}
#endif
__global__ void finish_seed(vdw_next::CompletionView v,State*s,const unsigned*base_error,u32 t_mode=0){if(threadIdx.x||blockIdx.x)return;
 if(*base_error){s->failed=1;return;}
 if constexpr(KIT_SINGLETON_SEED){if(s->coarse==v.r)s->best_c=s->best_f;}
 u32 c=1,f=1;
 if(!KIT_SINGLETON_SEED||s->coarse!=v.r) while(c+1<v.p){int z=v.at(c+1);if(z<0){s->failed=1;return;}if(z%(int)s->coarse)break;++c;}
 while(f+1<v.p){int z=v.at(f+1);if(z<0){s->failed=1;return;}if(z)break;++f;}
 if constexpr(KIT_SINGLETON_SEED){if(s->coarse==v.r)c=f;}
 s->prefix_c=c;s->prefix_f=f;
 u64 z=vdw_next::pack_best(c,1);if(z>s->best_c)s->best_c=z;z=vdw_next::pack_best(f,1);if(z>s->best_f)s->best_f=z;
 u32 q=s->best_f>>32;if(q>s->cap)q=s->cap;
 if(q<3){s->failed=2;return;}s->t=anchor_t(q,t_mode);
}
} // namespace kit
