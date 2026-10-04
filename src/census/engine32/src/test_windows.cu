#include <cuda_runtime.h>
#include "../vendor/research/baseline/baseline/previous/baseline/vendor/src/wheel_filter.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
using u32=std::uint32_t;
static void ck(cudaError_t e){if(e!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(e));}
template<class V>__global__ void windows(V v,int first,u32 count,u32*out){u32 j=blockIdx.x*blockDim.x+threadIdx.x;if(j>=count*4)return;u32 plane=j/count;int bit=first+(int)(j%count);out[j]=v.window(plane,bit)&v.valid(bit);}
static u32 label(u32 x){x^=x>>16;x*=0x7feb352du;x^=x>>15;x*=0x846ca68bu;return x^(x>>16);}
int main(){try{const u32 B=4096,np=4,flip=9;u32 checks=0;
 for(u32 G:{1u,2u,3u,6u}){u32 c=G==1?0:1,d=G==1?0:G-1,maximum=(B-c)/G,mirror_maximum=(B-d)/G,stride=std::max(maximum,mirror_maximum)/32+5;
 std::vector<u32>a(np*stride),b(np*stride);auto fill=[&](std::vector<u32>&out,u32 residue,u32 n){for(u32 x=0;x<=n;++x){u32 v=G==1&&x==0?0:label(G*x+residue);for(u32 p=0;p<np;++p)if((v>>p)&1)out[p*stride+1+x/32]|=1u<<(x%32);}};fill(a,c,maximum);fill(b,d,mirror_maximum);
 u32*da{},*db{},*dz{};ck(cudaMalloc(&da,a.size()*4));ck(cudaMalloc(&db,b.size()*4));ck(cudaMemcpy(da,a.data(),a.size()*4,cudaMemcpyHostToDevice));ck(cudaMemcpy(db,b.data(),b.size()*4,cudaMemcpyHostToDevice));
 int first=-(int)mirror_maximum-32;u32 count=maximum-first+33;ck(cudaMalloc(&dz,(size_t)np*count*4));
 if(G==1)windows<<<(np*count+127)/128,128>>>(vdw_next::View{da,stride,B,np,flip},first,count,dz);
 else windows<<<(np*count+127)/128,128>>>(vdw_next::ResidueView{da,db,stride,np,flip,maximum,mirror_maximum},first,count,dz);
 std::vector<u32>z((size_t)np*count);ck(cudaDeviceSynchronize());ck(cudaMemcpy(z.data(),dz,z.size()*4,cudaMemcpyDeviceToHost));
 for(u32 p=0;p<np;++p)for(u32 j=0;j<count;++j){int start=first+j;u32 expected=0;
  for(u32 k=0;k<32;++k){int x=start+k;bool valid=G==1?(x!=0&&x>=-(int)B&&x<=(int)B):(x>=-(int)mirror_maximum-1&&x<=(int)maximum);
   if(!valid)continue;u32 v;
   if(G==1)v=label((u32)std::abs(x));else v=label(x<0?G*(u32)(-x-1)+d:G*(u32)x+c);
   u32 bit=((v>>p)&1)^((x<0?(flip>>p):0)&1);if(bit)expected|=1u<<k;++checks;
  }
  if(z[p*count+j]!=expected)throw std::runtime_error("independent bit-window reconstruction mismatch");
 }
 cudaFree(da);cudaFree(db);cudaFree(dz);
 }
 std::cout<<"PASS "<<checks<<" independent bit checks across positive, negative, zero, mixed and guarded boundaries\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
