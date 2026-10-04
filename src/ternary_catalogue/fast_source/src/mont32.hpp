#pragma once
#include <cstdint>
#if defined(__CUDACC__)
#define MHD __host__ __device__
#else
#define MHD
#endif
namespace kit {
MHD inline std::uint32_t mont(std::uint32_t a,std::uint32_t b,std::uint32_t p,std::uint32_t ni){std::uint64_t t=(std::uint64_t)a*b,mn=(std::uint64_t)((std::uint32_t)t*ni)*p,z=t+mn;std::uint32_t q=z>>32;return (z<t||q>=p)?q-p:q;}
inline std::uint32_t ninverse(std::uint32_t p){std::uint32_t x=1;for(int i=0;i<5;++i)x*=2u-p*x;return 0u-x;}
}
#undef MHD
