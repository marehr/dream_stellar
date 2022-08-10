#include <cassert>
#include <iostream>
#include <limits>
#include <span>
#include <vector>

#include "simd_minimiser_window.hpp"

#include "diagnostics.hpp"
#include "compute_minimiser_blocks.hpp"
#include "compute_minimiser_blocks_simd.hpp"
#include "compute_minimiser_blocks_test.hpp"

using namespace stellar::test;

#include <bitset>

template <typename simd_t>
bool simd_equal(simd_t value, simd_t expected)
{
    size_t equal_counter = 0;
    for (size_t i = 0; i < simd_length<simd_t>; ++i)
        equal_counter += (value[i] == expected[i] ? 1 : 0);
    return equal_counter == simd_length<simd_t>;
};

void transpose_matrix_32x4x4_sse4(std::array<int32x4_t, 4> & matrix)
{
    __m128 tmp3, tmp2, tmp1, tmp0;
    tmp0 = _mm_unpacklo_ps((__m128)matrix[0], (__m128)matrix[1]);
    tmp2 = _mm_unpacklo_ps((__m128)matrix[2], (__m128)matrix[3]);
    tmp1 = _mm_unpackhi_ps((__m128)matrix[0], (__m128)matrix[1]);
    tmp3 = _mm_unpackhi_ps((__m128)matrix[2], (__m128)matrix[3]);
    matrix[0] = (int32x4_t)_mm_movelh_ps(tmp0, tmp2);
    matrix[1] = (int32x4_t)_mm_movehl_ps(tmp2, tmp0);
    matrix[2] = (int32x4_t)_mm_movelh_ps(tmp1, tmp3);
    matrix[3] = (int32x4_t)_mm_movehl_ps(tmp3, tmp1);
}

template <auto transpose_matrix_fn>
void transpose_matrix128_epi32_test()
{
    {
        std::array<int32x4_t, 4> matrix{};
        transpose_matrix_fn(matrix);
        assert(simd_equal(matrix[0], int32x4_t{}));
        assert(simd_equal(matrix[1], int32x4_t{}));
        assert(simd_equal(matrix[2], int32x4_t{}));
        assert(simd_equal(matrix[3], int32x4_t{}));
    }

    {
        std::array<int32x4_t, 4> matrix
        {
            int32x4_t{ 0,  1,  2,  3},
            int32x4_t{ 4,  5,  6,  7},
            int32x4_t{ 8,  9, 10, 11},
            int32x4_t{12, 13, 14, 15}
        };
        transpose_matrix_fn(matrix);
        assert(simd_equal(matrix[0], int32x4_t{0, 4, 8, 12}));
        assert(simd_equal(matrix[1], int32x4_t{1, 5, 9, 13}));
        assert(simd_equal(matrix[2], int32x4_t{2, 6, 10, 14}));
        assert(simd_equal(matrix[3], int32x4_t{3, 7, 11, 15}));
    }
}

void shuffle_test()
{
    __m256i a = _mm256_set_epi32(100, 101, 102, 103, 104, 105, 106, 107);
    __m256i idx = _mm256_set_epi32(4, 5, 6, 7, 0, 1, 2, 3);
    __m256i b = {};
    b = _mm256_permutevar8x32_epi32(a, idx);

    print_simd(b, 8);
}

static constexpr std::array<uint64_t, 256> _index_sequence_lookup{[]()
{
    std::array<uint64_t, 256> values{};

    auto bit_is_set = [](uint8_t const x, uint8_t const i) -> bool
    {
        return (x & (1 << i)) != 0;
    };

    for (uint16_t bit_sequence = 0; bit_sequence < 256; ++bit_sequence)
    {
        uint64_t index_sequence{0};

        for (uint8_t i = 8; i > 0;)
        {
            --i;
            if (bit_is_set(bit_sequence, i))
            {
                uint8_t new_value = i;
                index_sequence *= 256;
                index_sequence += new_value;
            }
        }
        values[bit_sequence] = index_sequence;
    }

    return values;
}()};

void _mask_compress_lookup_test()
{
    assert(_index_sequence_lookup[0b0000'0000u] == uint64_t{});

    assert(_index_sequence_lookup[0b0000'0001u] == uint64_t{0x00ull});
    assert(_index_sequence_lookup[0b0000'0010u] == uint64_t{0x01ull});
    assert(_index_sequence_lookup[0b0000'0100u] == uint64_t{0x02ull});
    assert(_index_sequence_lookup[0b0000'1000u] == uint64_t{0x03ull});
    assert(_index_sequence_lookup[0b0001'0000u] == uint64_t{0x04ull});
    assert(_index_sequence_lookup[0b0010'0000u] == uint64_t{0x05ull});
    assert(_index_sequence_lookup[0b0100'0000u] == uint64_t{0x06ull});
    assert(_index_sequence_lookup[0b1000'0000u] == uint64_t{0x07ull});

    // increasing sequences
    assert(_index_sequence_lookup[0b0000'0101u] == uint64_t{0x02'00ull});
    assert(_index_sequence_lookup[0b0001'0101u] == uint64_t{0x04'02'00ull});
    assert(_index_sequence_lookup[0b1001'0101u] == uint64_t{0x07'04'02'00ull});

    assert(_index_sequence_lookup[0b1110'1101u] == uint64_t{0x07'06'05'03'02'00ull});
    assert(_index_sequence_lookup[0b0000'0111u] == uint64_t{0x02'01'00ull});
}

int32x8_t simd_expand(/*uint8x8_t*/ uint64_t const src)
{
    __m128i src128 = _mm_cvtsi64_si128(src);
    int32x8_t dst{(int32x8_t)_mm256_cvtepu8_epi32(src128)};
    // same as
    // uint8x8_t src = (uint8x8_t)a;
    // int32x8_t dst = src;
    // for (int i = 0; i < 0; ++i)
    //     dst[i] = src[i];
    // _mm256_cvtepu8_epi32();
    return dst;
}

uint8_t simd_mask(int32x8_t const mask)
{
    uint32_t nearly_mask = _mm256_movemask_epi8((__m256i)mask);
    // std::cout << std::bitset<32>(nearly_mask) << std::endl;
    uint8_t final_mask{};

    // final_mask = (nearly_mask & 0x1) |
    //              (((nearly_mask & 0x10) >> 3) & (1 << 1)) |
    //              (((nearly_mask & 0x100) >> 6) & (1 << 2)) |
    //              (((nearly_mask & 0x1000) >> 9) & (1 << 3)) |
    //              (((nearly_mask & 0x10000) >> 12) & (1 << 4)) |
    //              (((nearly_mask & 0x100000) >> 15) & (1 << 5)) |
    //              (((nearly_mask & 0x1000000) >> 18) & (1 << 6)) |
    //              (((nearly_mask & 0x10000000) >> 21) & (1 << 7));
    // std::cout << std::bitset<32>(nearly_mask) << std::endl;
    // 0b1111'0000'1111'0000'1111'1111'1111'0000
    // 0b0001'0000'0001'0000'0001'0001'0001'0000
    uint32_t step0 = nearly_mask & 0b00010001'00010001'00010001'00010001;
    // std::cout << std::bitset<32>(step0) << std::endl;
    // 0b0000'0000'0100'0000'0100'0000'0100'0100
    uint32_t step1 = step0 >> 6;
    // std::cout << std::bitset<32>(step1) << std::endl;
    // 0b0001'0000'0001'0000'0001'0001'0001'0000
    // 0b0000'0000'0100'0000'0100'0000'0100'0100
    // | ---------------------------------------
    // 0b0001'0000'0101'0000'0101'0001'0101'0100
    uint32_t step2 = step0 | step1;
    // std::cout << std::bitset<32>(step2) << std::endl;
    // 0b0000'0010'0000'1010'0000'1010'0010'1010
    uint32_t step3 = step2 >> 3;
    // std::cout << std::bitset<32>(step3) << std::endl;
    // 0b0001'0000'0101'0000'0101'0001'0101'0100
    // 0b0000'0010'0000'1010'0000'1010'0010'1010
    // | ---------------------------------------
    // 0b0001'0010'0101'1010'0101'1011'0111'1110
    uint32_t step4 = step2 | step3;
    // std::cout << std::bitset<32>(step4) << std::endl;
    // 0b0001'0010'0101'1010'0101'1011'0111'1110
    // 0b0000'0000'0000'1111'0000'0000'0000'1111
    // & ---------------------------------------
    // 0b0000'0000'0000'1010'0000'0000'0000'1110
    uint32_t step5 = step4 & 0b00000000'00001111'00000000'00001111;
    // std::cout << std::bitset<32>(step5) << std::endl;
    // 0b0000'0000'0000'0000'0000'0000'1010'0000
    uint32_t step6 = step5 >> 12;
    // std::cout << std::bitset<32>(step6) << std::endl;
    // 0b0000'0000'0000'1010'0000'0000'0000'1110
    // 0b0000'0000'0000'0000'0000'0000'1010'0000
    // | ---------------------------------------
    // 0b0000'0000'0000'1010'0000'0000'1010'1110
    uint32_t step7 = step5 | step6;
    // std::cout << std::bitset<32>(step7) << std::endl;
    // 0b1010'1110
    uint8_t step8 = step7;
    // std::cout << std::bitset<32>(step8) << std::endl;
    final_mask = step8;

    // std::cout << std::bitset<32>(nearly_mask & 0x1) << std::endl;
    // std::cout << std::bitset<32>(nearly_mask & 0x10) << ", ";
    // std::cout << std::bitset<32>((nearly_mask & 0x10) >> 3) << ", ";
    // std::cout << std::bitset<32>(((nearly_mask & 0x10) >> 3) & 0b10) << std::endl;
    // std::cout << std::bitset<32>(nearly_mask & 0x100) << ", ";
    // std::cout << std::bitset<32>((nearly_mask & 0x100) >> 6) << ", ";
    // std::cout << std::bitset<32>(((nearly_mask & 0x100) >> 6) & 0b100) << std::endl;
    // std::cout << std::bitset<32>(nearly_mask & 0x1000) << ", ";
    // std::cout << std::bitset<32>((nearly_mask & 0x1000) >> 9) << ", ";
    // std::cout << std::bitset<32>(((nearly_mask & 0x1000) >> 9) & 0b1000) << std::endl;
    // std::cout << std::bitset<32>(nearly_mask & 0x10000) << ", ";
    // std::cout << std::bitset<32>((nearly_mask & 0x10000) >> 12) << ", ";
    // std::cout << std::bitset<32>(((nearly_mask & 0x10000) >> 12) & 0b10000) << std::endl;
    // std::cout << std::bitset<32>(nearly_mask & 0x100000) << ", ";
    // std::cout << std::bitset<32>((nearly_mask & 0x100000) >> 15) << ", ";
    // std::cout << std::bitset<32>(((nearly_mask & 0x100000) >> 15) & 0b100000) << std::endl;
    // for (int i = 8; i > 0;)
    // {
    //     --i;
    //     uint32_t bit_pos = 1 << (4 * i);
    //
    //     std::cout << "nearly_mask: " << std::bitset<32>(nearly_mask) << std::endl;
    //     std::cout << "bit_pos:    &" << std::bitset<32>(bit_pos) << std::endl;
    //     std::cout << "-------------" << std::bitset<32>(nearly_mask & bit_pos) << std::endl;
    //
    //     final_mask = final_mask << 1;
    //     final_mask |= ((nearly_mask & bit_pos) != 0 ? 1 : 0);
    //     std::cout << "final_mask:  " << std::bitset<32>(final_mask) << std::endl;
    //     // nearly_mask >= 1;
    // }
    // unsigned new_mask = nearly_mask / (unsigned)0b1111u;
    // std::cout << std::bitset<32>(final_mask) << std::endl;

    return final_mask;
}

void simd_mask_test()
{
    int32x8_t simd_four = int32x8_t{4, 4, 4, 4, 4, 4, 4, 4};
    assert(0b0000'0000 == simd_mask(int32x8_t{} == simd_four));
    assert(0b0000'0001 == simd_mask(int32x8_t{4} == simd_four));
    assert(0b0000'0010 == simd_mask(int32x8_t{0, 4} == simd_four));
    assert(0b0000'0100 == simd_mask(int32x8_t{0, 0, 4} == simd_four));
    assert(0b0000'1000 == simd_mask(int32x8_t{0, 0, 0, 4} == simd_four));
    assert(0b0001'0000 == simd_mask(int32x8_t{0, 0, 0, 0, 4} == simd_four));
    assert(0b0010'0000 == simd_mask(int32x8_t{0, 0, 0, 0, 0, 4} == simd_four));
    assert(0b0100'0000 == simd_mask(int32x8_t{0, 0, 0, 0, 0, 0, 4} == simd_four));
    assert(0b1000'0000 == simd_mask(int32x8_t{0, 0, 0, 0, 0, 0, 0, 4} == simd_four));
    assert(0b1111'1111 == simd_mask(int32x8_t{4, 4, 4, 4, 4, 4, 4, 4} == simd_four));
    assert(0b1010'1010 == simd_mask(int32x8_t{0, 4, 0, 4, 0, 4, 0, 4} == simd_four));
    assert(0b0101'0101 == simd_mask(int32x8_t{4, 0, 4, 0, 4, 0, 4, 0} == simd_four));
    assert(0b0101'0111 == simd_mask(int32x8_t{4, 4, 4, 0, 4, 0, 4, 0} == simd_four));
    assert(0b1010'1110 == simd_mask(int32x8_t{0, 4, 4, 4, 2, 4, 3, 4} == simd_four));
}

int32x8_t simd_mask_compress(uint8_t const mask, int32x8_t const src)
{
    // __m256i _mm256_maskz_compress_epi32(__mmask8 k, __m256i a) // swizzle
    // lookup index sequence for mask
    uint64_t idx_seq = _index_sequence_lookup[mask];
    int32x8_t idx_vec = simd_expand(idx_seq);

    return (int32x8_t)_mm256_permutevar8x32_epi32((__m256i)src, (__m256i)idx_vec);
}


// AVX512F + AVX512VL
void simd_mask_compress_test()
{
    int32x8_t values{100, 200, 300, 400, 500, 600, 700, 800};

    assert(simd_equal(simd_mask_compress(0b0000'0000, values), int32x8_t{100, 100, 100, 100, 100, 100, 100, 100}));
    assert(simd_equal(simd_mask_compress(0b0000'0001, values), int32x8_t{100, 100, 100, 100, 100, 100, 100, 100}));
    assert(simd_equal(simd_mask_compress(0b0000'0010, values), int32x8_t{200, 100, 100, 100, 100, 100, 100, 100}));
    assert(simd_equal(simd_mask_compress(0b0000'0100, values), int32x8_t{300, 100, 100, 100, 100, 100, 100, 100}));
    assert(simd_equal(simd_mask_compress(0b0000'1000, values), int32x8_t{400, 100, 100, 100, 100, 100, 100, 100}));
    assert(simd_equal(simd_mask_compress(0b0001'0000, values), int32x8_t{500, 100, 100, 100, 100, 100, 100, 100}));
    assert(simd_equal(simd_mask_compress(0b0010'0000, values), int32x8_t{600, 100, 100, 100, 100, 100, 100, 100}));
    assert(simd_equal(simd_mask_compress(0b0100'0000, values), int32x8_t{700, 100, 100, 100, 100, 100, 100, 100}));
    assert(simd_equal(simd_mask_compress(0b1000'0000, values), int32x8_t{800, 100, 100, 100, 100, 100, 100, 100}));
    assert(simd_equal(simd_mask_compress(0b0101'0101, values), int32x8_t{100, 300, 500, 700, 100, 100, 100, 100}));
    assert(simd_equal(simd_mask_compress(0b1010'1010, values), int32x8_t{200, 400, 600, 800, 100, 100, 100, 100}));
    assert(simd_equal(simd_mask_compress(0b1110'1010, values), int32x8_t{200, 400, 600, 700, 800, 100, 100, 100}));
    assert(simd_equal(simd_mask_compress(0b1110'1011, values), int32x8_t{100, 200, 400, 600, 700, 800, 100, 100}));
    // constexpr
}

int main()
{
    std::cout << "simd_mask_test" << std::endl;
    simd_mask_test();
    std::cout << "shuffle_test" << std::endl;
    shuffle_test();
    std::cout << "simd_mask_compress_test" << std::endl;
    simd_mask_compress_test();
    std::cout << "_mask_compress_lookup_test" << std::endl;
    _mask_compress_lookup_test();
    std::cout << "transpose_matrix128_epi32_test<transpose_matrix_32x4x4_sse4>" << std::endl;
    transpose_matrix128_epi32_test<transpose_matrix_32x4x4_sse4>();

    // gather_test();
    // compute_minimiser_blocks_test<compute_forward_full, compute_backward_full>();
    // // compute_minimiser_blocks_test<compute_forward_simd, compute_backward_simd>();
    // minimiser_test();
}
