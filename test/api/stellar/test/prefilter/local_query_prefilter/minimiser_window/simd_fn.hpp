
#pragma once

#include <x86intrin.h>

template <typename simd_t>
constexpr size_t simd_length = sizeof(simd_t) / sizeof(int);

constexpr size_t simd_len = 8;
using int32x8_t [[gnu::vector_size(sizeof(int) * simd_len)]] = int;
using int32x8_mask_t = int32x8_t;

using int32x4_t [[gnu::vector_size(sizeof(int) * 4)]] = int;
using int32x4_mask_t = int32x4_t;

template <typename simd_t>
inline simd_t simd_set(int offset)
{
    return simd_t{offset, offset, offset, offset, offset, offset, offset, offset};
}

inline int32x8_t simd_gather(size_t const window_size, size_t const offset, int const * __restrict source_ptr)
{
    int32x8_t r{};
    int32x8_t vindex{0, 1, 2, 3, 4, 5, 6, 7};
    vindex = vindex * (int)window_size + (int)offset;
    r = (int32x8_t)_mm256_i32gather_epi32(source_ptr, (__m256i)vindex, sizeof(int));
    return r;
}

template <typename simd_t>
simd_t simd_load(int const * mem_addr)
{
    struct missing_case{};

    if constexpr (std::is_same_v<simd_t, int32x8_t>)
    {
        return (simd_t)_mm256_loadu_si256(reinterpret_cast<__m256i_u const *>(mem_addr));
    } else
    {
        static_assert(std::is_same_v<simd_t, missing_case>);
        return {};
    }
}

inline void simd_store(int * mem_addr, int32x8_t a)
{
    _mm256_storeu_si256((__m256i *)mem_addr, (__m256i)a);
}

inline void simd_maskstore(int * mem_addr, int32x8_t mask,  int32x8_t a)
{
    _mm256_maskstore_epi32((int *)mem_addr, (__m256i)mask, (__m256i)a);
}

template <typename simd_t>
bool simd_equal(simd_t value, simd_t expected)
{
    size_t equal_counter = 0;
    for (size_t i = 0; i < simd_length<simd_t>; ++i)
        equal_counter += (value[i] == expected[i] ? 1 : 0);
    return equal_counter == simd_length<simd_t>;
};
