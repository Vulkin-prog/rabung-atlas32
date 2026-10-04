#pragma once
// Long-run catalogue capture for the singleton3 sparse scanner only.
namespace kit {
struct CapturedRun {u64 start;u32 length,colour;};
__device__ CapturedRun* capture_output;
__device__ unsigned long long* capture_count;
__device__ u32 capture_capacity;
__device__ u32 capture_minimum;
__global__ void capture_candidates(vdw_next::CompletionView v,State*s,const Candidate*in,unsigned long long capacity){
 v.t=s->t;const unsigned long long n=s->count;
 if(s->stop_flag||n>capacity||s->failed){if(!threadIdx.x&&!blockIdx.x)atomicExch(&s->failed,17u);return;}
 if(v.r!=3||s->coarse!=3||s->fine!=3){if(!threadIdx.x&&!blockIdx.x)atomicExch(&s->failed,18u);return;}
 u64 best=s->snap_c;
 for(unsigned long long i=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=(unsigned long long)gridDim.x*blockDim.x){
  auto word=in[i];u32 d=word.tagged_d&0x3fffffffu,tag=word.tagged_d>>30,g=tag==0?1:tag==1?2:tag==2?3:6;
  u32 bits=word.coarse|word.fine;
  while(bits){u32 bit=__ffs(bits)-1;bits&=bits-1;u64 wide=(u64)word.first_b+(u64)g*bit;
   if(!wide||wide>v.B||!d||d>v.D||d>v.B||d%g){atomicExch(&s->failed,19u);return;}
   u32 b=(u32)wide;if(!vdw_next::admissible(b,d,v.t))continue;
   u64 x=v.negative_quotient(b,d);if(!x||x+v.t>v.p){atomicExch(&s->failed,20u);return;}
   i64 a=-(i64)b;int colour=v.quotient_at(a,x,d);if(colour<0){atomicExch(&s->failed,21u);return;}
   u64 lo=x,hi=x+v.t-1;
   while(lo>1){int c=v.quotient_at(a+((i64)lo-1-(i64)x)*d,lo-1,d);if(c<0){atomicExch(&s->failed,22u);return;}if(c!=colour)break;--lo;}
   while(hi<v.p-1){int c=v.quotient_at(a+((i64)hi+1-(i64)x)*d,hi+1,d);if(c<0){atomicExch(&s->failed,23u);return;}if(c!=colour)break;++hi;}
   if(!packed_run(v,lo,hi,best)){atomicExch(&s->failed,24u);return;}
   u64 length=hi-lo+1;if(length>=capture_minimum){
    u64 mirror=v.p-hi,start=min(lo,mirror);unsigned long long j=atomicAdd(capture_count,1ull);
    if(j>=capture_capacity){atomicExch(&s->failed,25u);return;}
    // chi3(-1)=0, so folding preserves the observed colour.
    capture_output[j]={start,(u32)length,(u32)colour};
   }
  }
 }
 atomicMax(&s->best_c,(unsigned long long)best);atomicMax(&s->best_f,(unsigned long long)best);
}
} // namespace kit
