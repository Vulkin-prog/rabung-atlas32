#include "gpu_pipeline.cuh"
#include "small_table.cuh"
#include "seed_boundary.cuh"
static void ck(cudaError_t e){if(e!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(e));}
static u64 dense_best(const std::vector<unsigned char>&col,u32 p,u32 r,u32 B){
 u64 best=0;for(u32 lo=1;lo<p;){u32 hi=lo;while(hi+1<p&&col[hi+1]%r==col[lo]%r)++hi;
 if(lo<=B){u32 st=std::min(lo,p-hi);u64 z=vdw_next::pack_best(hi-lo+1,st);if(z>best)best=z;}lo=hi+1;}return best;
}
static void test(u32 p,u32 B,u32 c,u32 f,bool synthetic){
 std::vector<unsigned char>col(p);
 if(synthetic){u32 nh=(p-1)/2,sg=nh%f;for(u32 x=1;x<=nh;++x){u32 value=(x>=2040&&x<=2090)?0:(x-1)%2;col[x]=value;col[p-x]=(value+sg)%f;}}
 else{u32 g=root(p);u64 x=1;for(u32 i=0;i<p-1;++i){col[x]=i%f;x=x*g%p;}}
 unsigned char*dc{};kit::State*ds{},host{};unsigned*de{};ck(cudaMalloc(&dc,B+1));ck(cudaMalloc(&ds,sizeof(host)));ck(cudaMalloc(&de,4));
 ck(cudaMemcpy(dc,col.data(),B+1,cudaMemcpyHostToDevice));host.best_c=host.best_f=vdw_next::pack_best(1,1);host.cap=48;host.coarse=c;host.fine=f;host.D=1;host.t=3;
 ck(cudaMemcpy(ds,&host,sizeof(host),cudaMemcpyHostToDevice));ck(cudaMemset(de,0,4));
 vdw_next::CompletionView v{p,f,B,3,1,dc,nullptr};v.sigma_value=((p-1)/2)%f;
 u32 segments=(B+1023)/1024;kit::seed_runs<<<(segments+3)/4,128>>>(p,B,c,f,dc,ds);kit::seed_boundaries<<<((B-1)/1024+128)/128,128>>>(v,ds);kit::finish_seed<<<1,1>>>(v,ds,de);ck(cudaDeviceSynchronize());ck(cudaMemcpy(&host,ds,sizeof(host),cudaMemcpyDeviceToHost));require(!host.failed,"seed test failure");
 for(auto [r,packed]:std::vector<std::pair<u32,u64>>{{c,host.best_c},{f,host.best_f}}){u64 expected=dense_best(col,p,r,B),global=dense_best(col,p,r,p-1);
 require(packed>=expected&&packed<=global,"complete-seed inclusion/max bound");u32 len=packed>>32,lo=0xffffffffu-(u32)packed;require(lo&&((u64)lo+len<=p),"seed interval");
 for(u32 i=1;i<len;++i)require(col[lo+i]%r==col[lo]%r,"seed witness colour");if(lo>1)require(col[lo-1]%r!=col[lo]%r,"seed witness left boundary");if((u64)lo+len<p)require(col[lo+len]%r!=col[lo]%r,"seed witness right boundary");
 if(synthetic)require(packed==global,"warp-reduction regression (lane zero has no boundary run)");}
 cudaFree(dc);cudaFree(ds);cudaFree(de);
}
int main(){try{ck(cudaFree(nullptr));u32 count=0;
 for(auto [c,f]:std::vector<std::pair<u32,u32>>{{2,2},{2,4},{3,6},{5,5},{9,9}}){test(10081,10080,c,f,true);++count;}
 for(u32 p:{10009u,100003u,1000033u})for(auto [c,f]:kit::groups(p,"all"))for(u32 B:{1023u,1024u,1025u,2047u,2048u,2049u})if((u64)B*B>p&&B<p){test(p,B,c,f,false);++count;}
 std::cout<<"PASS "<<count<<" exact-seed boundary cases, independent dense colours and warp lane-zero regression\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
