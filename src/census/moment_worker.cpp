// The qualified point-counting algorithm is included unchanged. This adds only
// bounded batch transport and exact dense counts below its certified domain.
#define main qualified_mestre_main
#include "theorie/elliptic4_mestre.cpp"
#undef main
#include <sstream>
#include <string>
#include <cstdint>

void le(U x,unsigned bytes){for(unsigned i=0;i<bytes;++i)std::cout.put(char(x>>(8*i)));}
U small_A4(U p){
    if(p<=3)return 0;
    std::vector<unsigned char> square(p,0);
    for(U x=1;x<p;++x)square[x*x%p]=1;
    U n=0;
    for(U x=1;x+3<p;++x)n+=(square[x]==square[x+1]&&square[x]==square[x+2]&&square[x]==square[x+3]);
    return n;
}
int main(){try{
    std::ios::sync_with_stdio(false);
    std::cout<<"{\"type\":\"ready\",\"algorithm\":\"qualified_mestre_A4\",\"maximum_batch\":4096}\n"<<std::flush;
    std::string line;
    while(std::getline(std::cin,line)){
        if(line=="QUIT")return 0;
        std::istringstream input(line);U id,n;long long audit;
        if(!(input>>id>>audit>>n)||n>4096||audit< -1||audit>=static_cast<long long>(n))throw std::runtime_error("batch header");
        std::vector<U> primes(n);U previous=0;
        for(auto &p:primes){
            if(!(input>>p)||p<=previous||p>=(1ull<<32)||!exact_prime64(p))throw std::runtime_error("ordered prime batch");
            previous=p;
        }
        std::string excess;if(input>>excess)throw std::runtime_error("request suffix");
        auto begin=std::chrono::steady_clock::now();std::vector<std::uint32_t> values;values.reserve(n);
        U sum=0,small=0,max_attempt=0;std::string certificate;
        for(U i=0;i<n;++i){
            U p=primes[i],a;Card c;
            if(p<1000){a=small_A4(p);++small;}
            else{c=cardinality(p);a=quadruples(p,c.size);max_attempt=std::max(max_attempt,c.attempts);}
            if(a>0xffffffffull||(p>2&&a>(p-3)/2))throw std::runtime_error("A4 range");
            values.push_back(a);sum+=a;
            if(static_cast<long long>(i)==audit){
                if(p<1000)certificate="{\"status\":\"PASS\",\"p\":"+std::to_string(p)+",\"r\":2,\"quadruples\":"+std::to_string(a)+",\"certificate\":\"small_dense\"}";
                else{std::ostringstream text;auto old=std::cout.rdbuf(text.rdbuf());print(p,c,0.);std::cout.rdbuf(old);certificate=text.str();}
            }
        }
        U nanoseconds=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-begin).count();
        std::cout.write("A4BATCH1",8);le(id,8);le(n,4);le(certificate.size(),4);le(nanoseconds,8);le(sum,8);
        le(max_attempt,4);le(small,4);le(n?primes.front():0,8);le(n?primes.back():0,8);
        for(auto a:values)le(a,4);
        std::cout.write(certificate.data(),certificate.size());std::cout.flush();
    }
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
