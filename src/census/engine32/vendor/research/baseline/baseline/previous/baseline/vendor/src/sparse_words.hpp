#pragma once
#include <cstdint>
// Plain data contract; host decoder and CUDA emitter use exactly these fields.
namespace vdw_next {
struct alignas(16) SparseWord {
 std::uint32_t tagged_d; // low 30 bits d; high 2 bits encode g=1,2,3,6
 std::uint32_t first_b;  // candidate bit k means b=first_b+g*k
 std::uint32_t left;
 std::uint32_t right;
};
static_assert(sizeof(SparseWord)==16,"sparse output contract");
}
