// Independent exhaustive half-interval scan. Natural residues and corrected
// FP64 quotient replace the sparse scanner's Montgomery and anchor algorithm.
#include "src/wide_host.hpp"
#include <cuda_runtime.h>
#include <sstream>
static void ck(cudaError_t e){if(e!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(e));}
struct Summary {u32 len,first_colour,last_colour,prefix,suffix,best,start,first_x;};
__device__ __forceinline__ u64 fp_power(u64 x,u64 e,u64 p,double invp){
 u64 y=x;int bit=63-__clzll((unsigned long long)e);
 for(--bit;bit>=0;--bit){y=wide33::float_mod(y,y,p,invp);if((e>>bit)&1)y=wide33::float_mod(y,x,p,invp);}return y;
}
__global__ void dense_chunk(u64 p,u64 nh,u32 mask,u32 m,u64 e,double invp,const u64*roots,const u32*labels,
 u64 first_segment,u32 segments,u32 segment_size,Summary*out,u32*error){
 extern __shared__ u64 storage[];u64*rv=storage;u32*rl=(u32*)(storage+m);
 for(u32 j=threadIdx.x;j<m;j+=blockDim.x){rv[j]=roots[j];rl[j]=labels[j];}__syncthreads();
 u32 s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=segments)return;
 u64 first=1+(first_segment+s)*segment_size,last=first+segment_size-1;if(last>nh)last=nh;
 Summary z[8]{};u32 current[8]{};
 for(u32 i=0;i<8;++i)if(mask&(1u<<i)){z[i].len=(u32)(last-first+1);z[i].first_x=(u32)first;z[i].last_colour=~0u;}
 for(u64 x=first;x<=last;++x){
  u64 y=fp_power(x,e,p,invp);u32 lo=0,hi=m;while(lo<hi){u32 mid=(lo+hi)/2;if(rv[mid]<y)lo=mid+1;else hi=mid;}
  if(lo>=m||rv[lo]!=y){atomicExch(error,1u);return;}u32 label=rl[lo];
  #pragma unroll
  for(u32 i=0;i<8;++i)if(mask&(1u<<i)){u32 c=label%(i+2),run=c==z[i].last_colour?current[i]+1:1;current[i]=run;
   if(x==first)z[i].first_colour=c;
   if(c==z[i].first_colour&&z[i].prefix==x-first)++z[i].prefix;
   if(run>z[i].best){z[i].best=run;z[i].start=(u32)(x-run+1);}z[i].last_colour=c;
  }
 }
 for(u32 i=0;i<8;++i){z[i].suffix=current[i];out[(size_t)s*8+i]=z[i];}
}
static Summary merge(Summary a,Summary b){
 if(!a.len)return b;if(!b.len)return a;require((u64)a.first_x+a.len==b.first_x,"oracle segment ordering");
 u64 size=(u64)a.len+b.len;require(size<=UINT32_MAX,"oracle half length overflow");
 Summary c=a;c.len=(u32)size;c.last_colour=b.last_colour;
 if(b.best>c.best||(b.best==c.best&&b.start<c.start)){c.best=b.best;c.start=b.start;}
 if(a.last_colour==b.first_colour){u64 crossed=(u64)a.suffix+b.prefix,start=(u64)a.first_x+a.len-a.suffix;require(crossed<=UINT32_MAX&&start<=UINT32_MAX,"oracle crossing overflow");
  if(crossed>c.best||(crossed==c.best&&start<c.start)){c.best=(u32)crossed;c.start=(u32)start;}
  if(a.prefix==a.len)c.prefix=a.len+b.prefix;if(b.suffix==b.len)c.suffix=b.len+a.suffix;
 }else c.suffix=b.suffix;
 // For an equal boundary with a partial suffix in b, the suffix is b's suffix.
 if(a.last_colour==b.first_colour&&b.suffix!=b.len)c.suffix=b.suffix;
 return c;
}
int main(int argc,char**argv){try{
 u64 p=0,limit=0;u32 segment_size=1024,chunk_segments=65536;std::string scope="all";
 for(int i=1;i<argc;++i){require(i+1<argc,"oracle missing value");std::string k=argv[i],v=argv[++i];if(k=="--prime")p=parse64(v.c_str());else if(k=="--scope")scope=v;else if(k=="--limit")limit=parse64(v.c_str());else if(k=="--segment")segment_size=parse(v.c_str());else if(k=="--chunk-segments")chunk_segments=parse(v.c_str());else throw std::runtime_error("oracle unknown option "+k);}
 require(p>=1000&&p<(1ull<<33)&&wide33::prime(p),"oracle domain");require(segment_size>=32&&segment_size<=4096&&chunk_segments>=1&&chunk_segments<=262144,"oracle geometry");
 u32 mask=kit::requested_mask(p,scope),m=1;for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2)))m=std::lcm(m,r);
 u64 g=root(p),unit=3;while(std::gcd(unit,p-1)!=1)unit+=2;u64 alternate=wide33::power128(g,unit,p),z=wide33::power128(alternate,(p-1)/m,p),cur=1;
 std::vector<std::pair<u64,u32>>pairs(m);for(u32 j=0;j<m;++j){pairs[j]={cur,j};cur=(unsigned __int128)cur*z%p;}std::sort(pairs.begin(),pairs.end());
 std::vector<u64>rv(m);std::vector<u32>rl(m);for(u32 j=0;j<m;++j){rv[j]=pairs[j].first;rl[j]=pairs[j].second;}
 u64 nh=(p-1)/2,coverage=limit?std::min(nh,limit):nh,nsegments=(coverage+segment_size-1)/segment_size;bool complete=coverage==nh;
 u64 *dr{};u32*dl{},*de{};Summary*ds{};ck(cudaMalloc(&dr,m*8));ck(cudaMalloc(&dl,m*4));ck(cudaMalloc(&de,4));ck(cudaMalloc(&ds,(size_t)chunk_segments*8*sizeof(Summary)));
 ck(cudaMemcpy(dr,rv.data(),m*8,cudaMemcpyHostToDevice));ck(cudaMemcpy(dl,rl.data(),m*4,cudaMemcpyHostToDevice));ck(cudaMemset(de,0,4));
 Summary result[8]{};std::vector<Summary>h((size_t)chunk_segments*8);double gpu_ms=0;auto start=std::chrono::steady_clock::now();cudaEvent_t a,b;ck(cudaEventCreate(&a));ck(cudaEventCreate(&b));
 for(u64 first=0;first<nsegments;first+=chunk_segments){u32 count=(u32)std::min<u64>(chunk_segments,nsegments-first);
  ck(cudaEventRecord(a));dense_chunk<<<(count+127)/128,128,(size_t)m*12>>>(p,coverage,mask,m,(p-1)/m,1.0/(double)p,dr,dl,first,count,segment_size,ds,de);ck(cudaGetLastError());ck(cudaEventRecord(b));ck(cudaEventSynchronize(b));float dt;ck(cudaEventElapsedTime(&dt,a,b));gpu_ms+=dt;
  ck(cudaMemcpy(h.data(),ds,(size_t)count*8*sizeof(Summary),cudaMemcpyDeviceToHost));u32 err;ck(cudaMemcpy(&err,de,4,cudaMemcpyDeviceToHost));require(!err,"oracle natural character lookup failed");
  for(u32 s=0;s<count;++s)for(u32 i=0;i<8;++i)if(mask&(1u<<i))result[i]=merge(result[i],h[(size_t)s*8+i]);
 }
 u32 maximum[10]{},prefix[10]{},witness[10]{};
 for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2))){auto s=result[r-2];require(s.len==coverage&&s.first_x==1,"oracle incomplete coverage");maximum[r]=s.best;prefix[r]=s.prefix;witness[r]=s.start;
  if(complete&&nh%r==0){u64 middle=2ull*s.suffix,st=nh-s.suffix+1;require(middle<=UINT32_MAX&&st<=UINT32_MAX,"oracle full run overflow");
   if(middle>maximum[r]||(middle==maximum[r]&&st<witness[r])){maximum[r]=(u32)middle;witness[r]=(u32)st;}}
  if(complete)certify_witness33(p,r,maximum[r],prefix[r],witness[r]);
 }
 if(complete&&p<=2000000){u32 dm[10]{},dp[10]{},dw[10]{};kit::dense_host(p,mask,dm,dp,dw);for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2)))require(maximum[r]==dm[r]&&prefix[r]==dp[r]&&witness[r]==dw[r],"exhaustive GPU / dense CPU mismatch");}
 double total=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
 std::cout<<std::setprecision(12)<<"{\"type\":\"oracle\",\"status\":\"PASS\",\"p\":"<<p<<",\"scope\":\""<<scope<<"\",\"complete\":"<<(complete?"true":"false")<<",\"coverage\":"<<coverage<<",\"half_length\":"<<nh<<",\"arithmetic\":\"natural_FP64_corrected\",\"primitive_root\":"<<alternate<<",\"original_root\":"<<g<<",\"root_exponent\":"<<unit<<",\"gpu_ms\":"<<gpu_ms<<",\"total_ms\":"<<total<<",\"profiles\":[";
 bool comma=false;for(u32 r=2;r<=9;++r)if(mask&(1u<<(r-2))){if(comma)std::cout<<',';comma=true;u64 K=std::max<u64>((u64)maximum[r]+1,nh%r==0?2ull*prefix[r]+2:(u64)prefix[r]+2);
 std::cout<<"{\"r\":"<<r<<",\"max\":"<<maximum[r]<<",\"prefix\":"<<prefix[r]<<",\"K\":"<<K<<",\"witness\":"<<witness[r]<<'}';}
 std::cout<<"]}\n";cudaFree(dr);cudaFree(dl);cudaFree(de);cudaFree(ds);cudaEventDestroy(a);cudaEventDestroy(b);return 0;
 }catch(const std::exception&e){std::cerr<<"DENSE_ORACLE_FAILURE: "<<e.what()<<'\n';return 2;}}
