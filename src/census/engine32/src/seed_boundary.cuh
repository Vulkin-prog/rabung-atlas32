#pragma once
// Complete every run crossing a 1024-entry seed boundary or the table end.
// Together with per-segment maxima and the exact prefix, this gives the best
// full run having a t-anchor entirely inside [1,B]. Reflection covers d=1.
namespace kit {
__device__ __forceinline__ bool seed_pack(vdw_next::CompletionView v,u64 lo,u64 hi,u64&best){
 if(!lo||hi<lo||hi>=v.p)return false;u64 length=hi-lo+1,mirror=(u64)v.p-hi,start=lo<mirror?lo:mirror;
 if(length>0xffffffffull||!start||start>0xffffffffull)return false;
 u64 z=vdw_next::pack_best((u32)length,(u32)start);if(z>best)best=z;return true;
}
__global__ void seed_boundaries(vdw_next::CompletionView v,State*s){
 const u32 job=blockIdx.x*blockDim.x+threadIdx.x,count=(v.B-1)/1024+1;
 u64 bc=0,bf=0;bool ok=true;
 if(job<count){
  u64 x=job+1==count?v.B:(job+1)*1024ull;
  int left=v.at(x),right=x+1<v.p?v.at(x+1):-2;
  if(left<0||right==-1)ok=false;
  else if(right>=0){
   u32 coarse=s->coarse;int co=left>=(int)coarse?left-coarse:left,other=right>=(int)coarse?right-coarse:right;
   if(co==other){
    u64 lo=x,hi=x+1;
    while(lo>1){int z=v.at(lo-1);if(z<0){ok=false;break;}int c=z>=(int)coarse?z-coarse:z;if(c!=co)break;--lo;}
    while(ok&&hi+1<v.p){int z=v.at(hi+1);if(z<0){ok=false;break;}int c=z>=(int)coarse?z-coarse:z;if(c!=co)break;++hi;}
    if(ok)ok=seed_pack(v,lo,hi,bc);
    if(ok&&coarse==v.r)bf=bc;
    else if(ok){
     int previous=-1;u32 run=0;
     for(u64 pos=lo;pos<=hi;++pos){int z=v.at(pos);if(z<0){ok=false;break;}if(z==previous){if(run==0xffffffffu){ok=false;break;}++run;}else run=1;previous=z;
      if(!seed_pack(v,pos-run+1,pos,bf)){ok=false;break;}
     }
    }
   }
  }
 }
 if(!__all_sync(FULL,ok)){if(!(threadIdx.x&31))atomicExch(&s->failed,1u);return;}
 u32 hc=__reduce_max_sync(FULL,(u32)(bc>>32)),hf=__reduce_max_sync(FULL,(u32)(bf>>32));
 u32 lc=__reduce_max_sync(FULL,(u32)(bc>>32)==hc?(u32)bc:0u),lf=__reduce_max_sync(FULL,(u32)(bf>>32)==hf?(u32)bf:0u);
 if(!(threadIdx.x&31)){if(hc)atomicMax(&s->best_c,((unsigned long long)hc<<32)|lc);if(hf)atomicMax(&s->best_f,((unsigned long long)hf<<32)|lf);}
}
}
