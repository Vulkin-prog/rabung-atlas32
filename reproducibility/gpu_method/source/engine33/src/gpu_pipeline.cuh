#pragma once
#include <cuda_runtime.h>
#ifndef KIT_SHARE_PLANS
#define KIT_SHARE_PLANS 0
#endif
#ifndef KIT_COMPACT_FINE
#define KIT_COMPACT_FINE 0
#endif
#if KIT_COMPACT_FINE && (!KIT_ALIAS_READS || !KIT_WARP_PACK || KIT_FINE_SHARED)
#error Compact fine planes require parent aliasing and warp packing, with no full fine cache
#endif
#ifndef KIT_ALIAS_READS
#define KIT_ALIAS_READS 0
#endif
#if KIT_ALIAS_READS && (!KIT_STATIC_FINE_PLANES || !KIT_STATIC_RESIDUAL)
#error Parent aliasing needs the qualified fixed paired profiles
#endif
#ifndef KIT_SMEM_SPILL
#define KIT_SMEM_SPILL 0
#endif
#ifndef KIT_COPY_MODE
#define KIT_COPY_MODE 0
#endif
#if KIT_COPY_MODE==2
#include <cuda_pipeline.h>
#endif
#include "wide_host.hpp"
#include "../vendor/research/baseline/baseline/previous/src/warp_select.hpp"
#ifndef KIT_EXACT_TABLE_SEED
#define KIT_EXACT_TABLE_SEED 0
#endif
#ifndef KIT_SKIP_D1
#define KIT_SKIP_D1 0
#endif
static_assert(!KIT_SKIP_D1||KIT_EXACT_TABLE_SEED,"skip d=1 needs a complete table seed");
#ifndef KIT_SINGLETON_FULL_Q
#define KIT_SINGLETON_FULL_Q 0
#endif
#ifndef KIT_SINGLETON_FAST
#define KIT_SINGLETON_FAST 0
#endif
#ifndef KIT_STATIC_RESIDUAL
#define KIT_STATIC_RESIDUAL 0
#endif
#ifndef KIT_STATIC_FINE_PLANES
#define KIT_STATIC_FINE_PLANES 0
#endif
static_assert(!KIT_STATIC_RESIDUAL||KIT_STATIC_FINE_PLANES,"residual specialization requires fixed profile planes");

#ifndef KIT_FINE_SHARED
#define KIT_FINE_SHARED 0
#endif
#ifndef KIT_WARP_COMPLETION
#define KIT_WARP_COMPLETION 0
#endif
#ifndef KIT_CHUNK_SHIFT
#define KIT_CHUNK_SHIFT 0
#endif
#ifndef KIT_LOW_INTERLEAVE
#define KIT_LOW_INTERLEAVE 0
#endif
#ifndef KIT_LOW_CHUNK_SHIFT
#define KIT_LOW_CHUNK_SHIFT 5
#endif
#ifndef KIT_MID_CHUNK_SHIFT
#define KIT_MID_CHUNK_SHIFT 3
#endif
#ifndef KIT_HIGH_CHUNK_SHIFT
#define KIT_HIGH_CHUNK_SHIFT 2
#endif
#ifndef KIT_LOW_CYCLES
#define KIT_LOW_CYCLES 16
#endif
#ifndef KIT_MID_DEN_DIV
#define KIT_MID_DEN_DIV 2
#endif
static_assert(KIT_LOW_CHUNK_SHIFT>=0&&KIT_LOW_CHUNK_SHIFT<=7,"low chunk bound");
static_assert(KIT_MID_CHUNK_SHIFT>=0&&KIT_MID_CHUNK_SHIFT<=7,"mid chunk bound");
static_assert(KIT_HIGH_CHUNK_SHIFT>=0&&KIT_HIGH_CHUNK_SHIFT<=7,"high chunk bound");
static_assert(KIT_LOW_CYCLES>=1&&KIT_LOW_CYCLES<=128&&KIT_MID_DEN_DIV>=1&&KIT_MID_DEN_DIV<=16,"tier bounds");
#ifndef KIT_ADAPTIVE_CHUNKS
#define KIT_ADAPTIVE_CHUNKS 0
#endif
static_assert(KIT_CHUNK_SHIFT>=0&&KIT_CHUNK_SHIFT<=5,"bounded split");
namespace kit {
struct GInput{const u32*data;u32 stride,np,B,flip;};
constexpr u32 FULL=0xffffffffu;
struct WP{vdw_periodic::DevicePlan p;u32 end;};
// Plain storage has no implicit initialization in CUDA shared memory.
struct SharedGeom{u32 bf,bl,ef,el,F,L,step,middle,bands,P,cycle,mid,total;float inverse;u32 end[32];};
struct SharedPlan{SharedGeom a,b;};
constexpr size_t PLAN_CACHE_BYTES=KIT_SHARE_PLANS?4*sizeof(SharedPlan):0;
static_assert(!KIT_SHARE_PLANS||(KIT_ADAPTIVE_CHUNKS&&KIT_LOW_CHUNK_SHIFT==5&&KIT_MID_CHUNK_SHIFT==3&&KIT_HIGH_CHUNK_SHIFT==3),"shared plans require the qualified 32/8/8 partition");
__device__ __forceinline__ void store_plan(SharedGeom &s,WP v,u32 lane){
 if(!lane){auto p=v.p;auto a=p.shape;s.bf=a.begin.first;s.bl=a.begin.last;s.ef=a.end.first;s.el=a.end.last;
 s.F=a.F;s.L=a.L;s.step=a.step;s.middle=a.middle;s.bands=a.bands;s.P=p.P;s.cycle=p.cycle;s.mid=p.mid_words;s.total=p.total;s.inverse=p.inverse_cycle;}
 s.end[lane]=v.end;
}
__device__ __forceinline__ WP load_plan(const SharedGeom&s,u32 lane){
 vdw_periodic::Shape a{{s.bf,s.bl},{s.ef,s.el},s.F,s.L,s.step,s.middle,s.bands};
 return {{a,s.P,s.cycle,s.mid,s.total,s.inverse},s.end[lane]};
}

__device__ __forceinline__ WP build_plan(vdw_periodic::Shape shape){
 u32 lane=threadIdx.x&31,P=shape.middle?vdw_periodic::period(shape.step):1;
 u32 sum=lane<P?(shape.middle?vdw_periodic::band_words(shape,lane):1):0;
 #pragma unroll
 for(u32 distance=1;distance<32;distance*=2){u32 before=__shfl_up_sync(FULL,sum,distance);if(lane>=distance)sum+=before;}
 u32 cycle=__shfl_sync(FULL,sum,P-1),tail=shape.middle%P;
 u32 tw=__shfl_sync(FULL,sum,tail?tail-1:0);if(!tail)tw=0;
 u32 mid=(shape.middle/P)*cycle+tw;
 float inverse=0;
#if KIT_RECIP_SELECT
 // Exact corrected quotient for integers <=2^22; otherwise preserve integer division.
 if(cycle<=(1u<<22)&&mid<=(1u<<22))inverse=__frcp_rn((float)cycle);
#endif
 return {{shape,P,cycle,mid,vdw_periodic::words(shape.begin)+mid+vdw_periodic::words(shape.end),inverse},sum};
}
struct Fetch{u32 value;__device__ __forceinline__ u32 operator()(u32 lane)const{return __shfl_sync(FULL,value,lane);}};
__device__ __forceinline__ vdw_periodic::Task select(WP a,WP b,u32 rank){
 u32 mask=__ballot_sync(FULL,rank<a.p.total);
 if(mask==FULL)return vdw_periodic::warp_select(a.p,rank,Fetch{a.end});
 if(!mask)return vdw_periodic::warp_select(b.p,rank-a.p.total,Fetch{b.end});
 auto x=vdw_periodic::warp_select(a.p,rank,Fetch{a.end});
 auto y=vdw_periodic::warp_select(b.p,rank<a.p.total?b.p.total:rank-a.p.total,Fetch{b.end});
 return rank<a.p.total?x:y;
}
template<int G> __device__ __forceinline__ auto row_view(GInput v,u32 row){
 if constexpr(G==1)return vdw_next::View{v.data,v.stride,v.B,v.np,v.flip};
 else{constexpr u32 N=G==2?1:2;u32 c=row==0?1:G-1;
  return vdw_next::ResidueView{v.data+row*v.np*v.stride,v.data+(N-1-row)*v.np*v.stride,v.stride,v.np,v.flip,(v.B-c)/G,(v.B-(G-c))/G};}
}
// Fine encodings share every non-residual plane with their coarse parent.
// 2->4: fine bit 0 = coarse bit 0; 3->6: fine bits 1,2 = coarse bits 0,1.
constexpr bool parent_plane_identity(){
 for(u32 c=0;c<4;++c)if((c&1)!=((c%2)&1))return false;
 for(u32 c=0;c<6;++c){u32 fine=(c&1)+2*(c%3);if((fine>>1)!=(c%3))return false;}return true;
}
static_assert(parent_plane_identity(),"fine-parent encoding identity");
template<class C,class F,int NP,int NF,int RES>struct ParentFine {
 C coarse;F fine;static constexpr u32 np=NF,residual=RES;
 __host__ __device__ __forceinline__ u32 parent_plane(u32 k)const{return RES==0?k-1:k;}
 __host__ __device__ __forceinline__ u32 aligned(u32 k,u32 w)const{return k==residual?fine.aligned(KIT_COMPACT_FINE?0:k,w):coarse.aligned(parent_plane(k),w);}
 __host__ __device__ __forceinline__ u32 positive(u32 k,int b)const{return k==residual?fine.positive(KIT_COMPACT_FINE?0:k,b):coarse.positive(parent_plane(k),b);}
 __host__ __device__ __forceinline__ u32 window(u32 k,int b)const{return k==residual?fine.window(KIT_COMPACT_FINE?0:k,b):coarse.window(parent_plane(k),b);}
 __host__ __device__ __forceinline__ u32 valid(int b)const{return coarse.valid(b);}
};
template<int NP,int NF,int RES,class C,class F>__device__ __forceinline__ auto parent_fine(C c,F f){
 if constexpr((NP==1&&NF==2)||(NP==2&&NF==3))return ParentFine<C,F,NP,NF,RES>{c,f};
 else return f;
}
__device__ __forceinline__ u32 anchor_t(u32 q,u32 mode){
 u32 half=(q+1)/2,t=mode==1?(q+2)/3:mode==2?(q+3)/4:mode==3?(2*q+2)/3:mode==4?(3*q+3)/4:mode==5?q:half;
 if(mode>=6&&mode<=9)t=mode==6?half-1:mode==7?half+1:mode==8?half-2:half+2;
 if(mode>=10){u32 space=q-half+1,lower=1u<<(31-__clz(space)),upper=lower==space?lower:2*lower;
  u32 chosen=mode==11?upper:mode==12?(space-lower<=upper-space?lower:upper):lower;
  u32 proposed=q>=chosen?q-chosen+1:3;
  if(mode!=13||(proposed>=half&&proposed-half<=2))t=proposed;
 }
 return min(q,max(3u,t));
}
__global__ void epoch_snapshot(State*s,u32 lo,u32 hi,bool dynamic_t=false,u32 t_mode=0){if(!threadIdx.x&&!blockIdx.x){
  if(dynamic_t){u32 q=min((u32)(s->best_c>>32),(u32)(s->best_f>>32));q=min(q,s->cap);
   u32 next=anchor_t(q,t_mode);s->t=max(s->t,next);}
  take_snapshot(*s,lo,hi);
 }}
// Each launch is bounded by one denominator epoch. Thresholds are immutable until
// all four groups and completion finish on the SAME stream/graph dependency.
template<int NP,int G,bool SHARED,int NF=0,int RES=(NP==1?1:0)>
__global__ void generate(GInput cp,GInput fp,State *state,u32 selected,bool target,
                         Candidate*out,unsigned long long capacity,bool stage_fine){
#if KIT_SMEM_SPILL
 // The L2 specialization has no dynamic-shared accesses after if-constexpr
 // elimination. The pragma is never used by a table-caching specialization.
 if constexpr(!SHARED)asm volatile(".pragma \"enable_smem_spilling\";");
#endif
 constexpr int NR=(G==1||G==2)?1:2;constexpr u32 TAG=G==1?0:G==2?1:G==3?2:3;
 extern __shared__ __align__(16) u32 smem[];
 if constexpr(SHARED){
#if KIT_COPY_MODE==0
  for(u32 j=threadIdx.x;j<NR*NP*cp.stride;j+=blockDim.x)smem[j]=cp.data[j];
#else
  const u32 total=NR*NP*cp.stride;
  for(u32 j=threadIdx.x*4;j<total;j+=blockDim.x*4){
   if(j+3<total){
#if KIT_COPY_MODE==1
    reinterpret_cast<int4*>(smem)[j/4]=reinterpret_cast<const int4*>(cp.data)[j/4];
#else
    __pipeline_memcpy_async(smem+j,cp.data+j,16);
#endif
   }else{
    for(u32 k=j;k<total;++k)smem[k]=cp.data[k];
   }
  }
#if KIT_COPY_MODE==2
  __pipeline_commit();__pipeline_wait_prior(0);
#endif
#endif
#if KIT_FINE_SHARED
  if constexpr(NP<=2){
   if(stage_fine){
    u32 offset=(NR*NP*cp.stride+3)&~3u;u32*dst=smem+offset;
    const u32 total=NR*fp.np*fp.stride;
    for(u32 j=threadIdx.x*4;j<total;j+=blockDim.x*4){
     if(j+3<total){
#if KIT_COPY_MODE==2
      __pipeline_memcpy_async(dst+j,fp.data+j,16);
#else
      reinterpret_cast<int4*>(dst)[j/4]=reinterpret_cast<const int4*>(fp.data)[j/4];
#endif
     }else for(u32 k=j;k<total;++k)dst[k]=fp.data[k];
    }
#if KIT_COPY_MODE==2
    __pipeline_commit();__pipeline_wait_prior(0);
#endif
    fp.data=dst;
   }
  }
#endif
  __syncthreads();cp.data=smem;
 }
 cp.np=NP;if constexpr(NF>0)fp.np=NF;
#if KIT_STATIC_RESIDUAL
 // The qualified paired groups are 2->4 (residual bit 1) and 3->6 (bit 0).
 if constexpr(NP==1&&NF==2)selected=2;
 else if constexpr(NP==2&&NF==3)selected=1u<<RES;
#endif
const u32 lane=threadIdx.x&31,warp=(blockIdx.x*blockDim.x+threadIdx.x)/32;
 const u32 nw=gridDim.x*(blockDim.x/32),lo=state->lo,hi=state->hi,t=state->t,q=state->q;
 const u32 qc=state->qc,qf=state->qf,space=q-t+1;
 const u32 cycle0=(lo-1)/6,cycles=(hi-1+5)/6-cycle0;
 const u32 per=(G==1||G==2)?2:1,original_slots=cycles*per*NR;
#if KIT_ADAPTIVE_CHUNKS
 // Exact denominator tiers and chunk counts; defaults: <96:32, lower half:8, upper half:4.
 // Piecewise integer rank mapping partitions each original word exactly once.
 const u32 j0=cycle0<KIT_LOW_CYCLES?min(original_slots,(KIT_LOW_CYCLES-cycle0)*per*NR):0;
 const u32 half_cycle=(state->D+6*KIT_MID_DEN_DIV-1)/(6*KIT_MID_DEN_DIV);
 const u32 j1=max(j0,cycle0<half_cycle?min(original_slots,(half_cycle-cycle0)*per*NR):0);
 const u32 edge0=j0<<KIT_LOW_CHUNK_SHIFT,edge1=edge0+((j1-j0)<<KIT_MID_CHUNK_SHIFT);
 const u32 slots=edge1+((original_slots-j1)<<KIT_HIGH_CHUNK_SHIFT);
#else
 const u32 slots=original_slots<<KIT_CHUNK_SHIFT;
#endif
#if KIT_SHARE_PLANS
 __shared__ SharedPlan shared_plans[4];
 const u32 local_warp=threadIdx.x/32;
 // Every physical thread visits the same batches, including inactive warps.
 for(u32 batch=blockIdx.x*(blockDim.x/32);batch<slots;batch+=nw){
  const u32 job=batch+local_warp;bool plan_active=job<slots;
  __syncthreads(); // previous readers have copied the plans to registers
#else
 for(u32 job=warp;job<slots;job+=nw){
#endif
  u32 original_job,chunk,chunk_shift;
#if KIT_ADAPTIVE_CHUNKS
  if(job<edge0){
#if KIT_LOW_INTERLEAVE
   // Spread the 32 chunks of each small denominator over separate CTAs.
   original_job=job%j0;chunk=job/j0;chunk_shift=KIT_LOW_CHUNK_SHIFT;
#else
   original_job=job>>KIT_LOW_CHUNK_SHIFT;chunk=job&((1u<<KIT_LOW_CHUNK_SHIFT)-1);chunk_shift=KIT_LOW_CHUNK_SHIFT;
#endif
  }
  else if(job<edge1){u32 rank=job-edge0;original_job=j0+(rank>>KIT_MID_CHUNK_SHIFT);chunk=rank&((1u<<KIT_MID_CHUNK_SHIFT)-1);chunk_shift=KIT_MID_CHUNK_SHIFT;}
  else{u32 rank=job-edge1;original_job=j1+(rank>>KIT_HIGH_CHUNK_SHIFT);chunk=rank&((1u<<KIT_HIGH_CHUNK_SHIFT)-1);chunk_shift=KIT_HIGH_CHUNK_SHIFT;}
#else
  original_job=job>>KIT_CHUNK_SHIFT;chunk=job&((1u<<KIT_CHUNK_SHIFT)-1u);chunk_shift=KIT_CHUNK_SHIFT;
#endif
  u32 row=original_job%NR,n=original_job/NR,d;
  if constexpr(G==1)d=6*(cycle0+n/2)+(n%2?5:1);
  if constexpr(G==2)d=6*(cycle0+n/2)+(n%2?4:2);
  if constexpr(G==3)d=6*(cycle0+n)+3;
  if constexpr(G==6)d=6*(cycle0+n)+6;
  u32 c=G==1?0:row==0?1:G-1,phase=((t-1)/2)%space;
#if KIT_SHARE_PLANS
  plan_active=plan_active&&d>=lo&&d<hi;
#if KIT_SKIP_D1
  if constexpr(G==1)plan_active=plan_active&&d!=1;
#endif
  // Low tier: one denominator per CTA; other tiers: four groups of 8 warps.
  const u32 cache_slot=chunk_shift==5?0:local_warp/8;
  if(plan_active&&(local_warp&((1u<<chunk_shift)-1))==0){
   auto shapes=vdw_periodic::make_shapes(cp.B,d,G,c,t,space,phase);
   auto pa=build_plan(shapes.a),pb=build_plan(shapes.b);
   store_plan(shared_plans[cache_slot].a,pa,lane);store_plan(shared_plans[cache_slot].b,pb,lane);
  }
  __syncthreads();
  if(!plan_active)continue;
  auto a=load_plan(shared_plans[cache_slot].a,lane),b=load_plan(shared_plans[cache_slot].b,lane);
#else
  if(d<lo||d>=hi)continue;
#if KIT_SKIP_D1
  if constexpr(G==1)if(d==1)continue;
#endif
  auto shapes=vdw_periodic::make_shapes(cp.B,d,G,c,t,space,phase);
  auto a=build_plan(shapes.a),b=build_plan(shapes.b);
#endif
  auto cv=row_view<G>(cp,row);
#if KIT_COMPACT_FINE
  GInput residual_input=fp;
  if constexpr((NP==1&&NF==2)||(NP==2&&NF==3)){residual_input.np=1;residual_input.flip=(fp.flip>>RES)&1;}
  auto raw_fv=row_view<G>(residual_input,row);
#else
  auto raw_fv=row_view<G>(fp,row);
#endif
#if KIT_ALIAS_READS
  auto fv=parent_fine<NP,NF,RES>(cv,raw_fv);
#else
  auto fv=raw_fv;
#endif
  u32 total=a.p.total+b.p.total;
  for(u32 start=chunk*32;start<total;start+=(32u<<chunk_shift)){auto task=select(a,b,start+lane);vdw_target::Masks masks{0,0};
   if(task.valid){u32 bits=vdw_next::filter_range(cv,d/G,t,task.word,task.first,task.last,nullptr,vdw_next::primitive_mask<G>(d,c,task.word));
    #if KIT_SINGLETON_FAST
    // For a valid divisor/refinement pair, equal nonzero plane counts imply
    // coarse==fine. Its unused refinement path can be removed at compile time.
    if constexpr(NF>0 && NF==NP){
     // Geometry still samples anchors for the certified capped q-window.
     // Each run as long as qc contains such a window. Reject only when both
     // observed mismatches certify that this anchor lies in a run shorter
     // than qc. Values outside the table stay unknown and admissible.
     const u32 screen_q=KIT_SINGLETON_FULL_Q?min(qc,256u):q;
     masks={vdw_cut::screen(cv,d/G,t,screen_q,task.word,bits,1,nullptr),0};
    }
    else
#endif
    {if(target&&selected)masks=vdw_target::split(cv,fv,selected,d/G,t,qc,qf,task.word,bits);
    else masks={vdw_cut::screen(cv,d/G,t,q,task.word,bits,1,nullptr),0};}

   }
   u32 votes=__ballot_sync(FULL,(masks.coarse|masks.fine)!=0);
   if(votes){unsigned long long at=0;if(!lane)at=atomicAdd(&state->count,(unsigned long long)__popc(votes));at=__shfl_sync(FULL,at,0);
    if(masks.coarse|masks.fine){u32 rank=__popc(votes&((1u<<lane)-1u));
     if(at+rank>=capacity)atomicExch(&state->failed,1u);
     else out[at+rank]={d|(TAG<<30),G*32*task.word+c,masks.coarse,masks.fine};}
   }
  }
 }
}
__global__ void complete_candidates(vdw_next::CompletionView v,State*s,const Candidate*in,unsigned long long capacity){
 v.t=s->t; const unsigned long long n=s->count;if(n>capacity||s->failed){if(!threadIdx.x&&!blockIdx.x)atomicExch(&s->failed,1u);return;}
 u64 bc=s->snap_c,bf=s->snap_f;bool ok=true;
 for(unsigned long long i=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=(unsigned long long)gridDim.x*blockDim.x)
  if(!kit::complete(v,s->coarse,in[i],bc,bf)){ok=false;break;}
#if KIT_WARP_COMPLETION
 // All physical lanes converge here, including lanes with no descriptor.
 if(!__all_sync(FULL,ok)){if(!(threadIdx.x&31))atomicExch(&s->failed,1u);return;}
 u32 hc=__reduce_max_sync(FULL,(u32)(bc>>32));
 u32 lc=__reduce_max_sync(FULL,(u32)(bc>>32)==hc?(u32)bc:0u);
 u32 hf=__reduce_max_sync(FULL,(u32)(bf>>32));
 u32 lf=__reduce_max_sync(FULL,(u32)(bf>>32)==hf?(u32)bf:0u);
 bc=((u64)hc<<32)|lc;bf=((u64)hf<<32)|lf;
 if(!(threadIdx.x&31)){
  if(bc>s->snap_c)atomicMax(&s->best_c,(unsigned long long)bc);
  if(bf>s->snap_f)atomicMax(&s->best_f,(unsigned long long)bf);
 }
#else
 if(!ok){atomicExch(&s->failed,1u);return;}
 atomicMax(&s->best_c,(unsigned long long)bc);atomicMax(&s->best_f,(unsigned long long)bf);
#endif
}
__global__ void end_epoch(State*s,unsigned long long capacity){if(!threadIdx.x&&!blockIdx.x){if(s->count>capacity)s->failed=1;else ++s->epochs_done;}}
struct Resources{int registers=0,local_bytes=0,active_blocks=0;size_t shared_bytes=0;bool used_shared=false,fine_cached=false;};
template<int NP,int G,bool S,int NF=0,int RES=(NP==1?1:0)> static cudaError_t launch_typed(GInput c,GInput f,State*s,u32 sel,bool target,Candidate*out,
 unsigned long long cap,int blocks,int threads,cudaStream_t stream,Resources*res=nullptr){
 if(KIT_SHARE_PLANS&&threads!=128&&threads!=256&&threads!=512&&threads!=1024)return cudaErrorInvalidValue;
 constexpr int N=G==1||G==2?1:2;size_t bytes=S?N*NP*(size_t)c.stride*4:0;bool stage_fine=false;
#if KIT_FINE_SHARED
 if constexpr(S && NP<=2){
  static const int limit=[](){int n=0;auto e=cudaDeviceGetAttribute(&n,cudaDevAttrMaxSharedMemoryPerBlockOptin,0);return e==cudaSuccess?n:0;}();
  size_t combined=((bytes+15)&~size_t(15))+N*(size_t)f.np*f.stride*4;
  stage_fine=target&&sel&&combined<=(size_t)limit;
  if(stage_fine)bytes=combined;
 }
#endif
 if(res){if constexpr(S){if(bytes+PLAN_CACHE_BYTES>48*1024){static const cudaError_t configured=[](){int limit=0;auto e=cudaDeviceGetAttribute(&limit,cudaDevAttrMaxSharedMemoryPerBlockOptin,0);if(e!=cudaSuccess)return e;
        cudaFuncAttributes a{};e=cudaFuncGetAttributes(&a,generate<NP,G,S,NF,RES>);if(e!=cudaSuccess)return e;
        return cudaFuncSetAttribute(generate<NP,G,S,NF,RES>,cudaFuncAttributeMaxDynamicSharedMemorySize,limit-(int)a.sharedSizeBytes);}();if(configured!=cudaSuccess)return configured;}}
 cudaFuncAttributes a{};auto e=cudaFuncGetAttributes(&a,generate<NP,G,S,NF,RES>);if(e!=cudaSuccess)return e;
  int active=0;e=cudaOccupancyMaxActiveBlocksPerMultiprocessor(&active,generate<NP,G,S,NF,RES>,threads,bytes);if(e!=cudaSuccess)return e;
  *res={a.numRegs,(int)a.localSizeBytes,active,bytes+PLAN_CACHE_BYTES,S,stage_fine};if(!active)return cudaErrorInvalidConfiguration;return cudaSuccess;}
 generate<NP,G,S,NF,RES><<<blocks,threads,bytes,stream>>>(c,f,s,sel,target,out,cap,stage_fine);return cudaPeekAtLastError();
}
template<int NP,int G,int NF=0,int RES=(NP==1?1:0)>static cudaError_t launch_choice(GInput c,GInput f,State*s,u32 sel,bool target,Candidate*out,unsigned long long cap,
 int blocks,int threads,bool shared,cudaStream_t stream,Resources*r){
 return shared?launch_typed<NP,G,true,NF,RES>(c,f,s,sel,target,out,cap,blocks,threads,stream,r):launch_typed<NP,G,false,NF,RES>(c,f,s,sel,target,out,cap,blocks,threads,stream,r);
}
static cudaError_t launch_group(int group,GInput c,GInput f,State*s,u32 sel,bool target,Candidate*out,unsigned long long cap,
 int blocks,int threads,bool shared,cudaStream_t stream,Resources*r=nullptr){
#if KIT_STATIC_RESIDUAL
 if((c.np==1&&f.np==2&&sel!=2)||(c.np==2&&f.np==3&&sel!=1&&sel!=4))return cudaErrorInvalidValue;
#endif
#if KIT_STATIC_FINE_PLANES
 // Only the current benchmark groups are instantiated: 2->2/4, 3->6,
 // and singleton 5/7/8/9. Unsupported combinations are refused explicitly.
 if(c.np==2&&f.np==3&&sel==4){
  switch(group){
   case 1:return launch_choice<2,1,3,2>(c,f,s,sel,target,out,cap,blocks,threads,shared,stream,r);
   case 2:return launch_choice<2,2,3,2>(c,f,s,sel,target,out,cap,blocks,threads,shared,stream,r);
   case 3:return launch_choice<2,3,3,2>(c,f,s,sel,target,out,cap,blocks,threads,shared,stream,r);
   case 6:return launch_choice<2,6,3,2>(c,f,s,sel,target,out,cap,blocks,threads,shared,stream,r);
  }
 }
 #define GO(N,G,F) return launch_choice<N,G,F>(c,f,s,sel,target,out,cap,blocks,threads,shared,stream,r)
 #define GROUP(N,F) switch(group){case 1:GO(N,1,F);case 2:GO(N,2,F);case 3:GO(N,3,F);case 6:GO(N,6,F);}
 switch(c.np){
  case 1:if(f.np==1){GROUP(1,1)}else if(f.np==2){GROUP(1,2)}break;
  case 2:if(f.np==2){GROUP(2,2)}else if(f.np==3){GROUP(2,3)}break;
  case 3:if(f.np==3){GROUP(3,3)}break;
  case 4:if(f.np==4){GROUP(4,4)}break;
 }
 #undef GROUP
 #undef GO
#else
 #define GO(N,G) return launch_choice<N,G>(c,f,s,sel,target,out,cap,blocks,threads,shared,stream,r)
 #define GROUP(N) switch(group){case 1:GO(N,1);case 2:GO(N,2);case 3:GO(N,3);case 6:GO(N,6);}
 switch(c.np){case 1:GROUP(1);case 2:GROUP(2);case 3:GROUP(3);case 4:GROUP(4);}
 #undef GROUP
 #undef GO
#endif
 return cudaErrorInvalidValue;
}
} // namespace kit
