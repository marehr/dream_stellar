
#pragma once

#include "../simd_swizzle_state_rules.hpp"

simd_swizzle_state_rules shuffle_epi32_transitions{
    // ------ shuffle_epi32 ------
    simd_swizzle_op{{0, 0, 0, 0,/*|*/ 4, 4, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 0, 0,/*|*/ 5, 4, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 0, 0,/*|*/ 6, 4, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 0, 0,/*|*/ 7, 4, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 0, 0,/*|*/ 4, 5, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 0, 0,/*|*/ 5, 5, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 0, 0,/*|*/ 6, 5, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 0, 0,/*|*/ 7, 5, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 0, 0,/*|*/ 4, 6, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 0, 0,/*|*/ 5, 6, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 0, 0,/*|*/ 6, 6, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 0, 0,/*|*/ 7, 6, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 0, 0,/*|*/ 4, 7, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 0, 0,/*|*/ 5, 7, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 0, 0,/*|*/ 6, 7, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 0, 0,/*|*/ 7, 7, 4, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0000'1111}, argument_none{}},
    simd_swizzle_op{{0, 0, 1, 0,/*|*/ 4, 4, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 1, 0,/*|*/ 5, 4, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 1, 0,/*|*/ 6, 4, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 1, 0,/*|*/ 7, 4, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 1, 0,/*|*/ 4, 5, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 1, 0,/*|*/ 5, 5, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 1, 0,/*|*/ 6, 5, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 1, 0,/*|*/ 7, 5, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 1, 0,/*|*/ 4, 6, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 1, 0,/*|*/ 5, 6, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 1, 0,/*|*/ 6, 6, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 1, 0,/*|*/ 7, 6, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 1, 0,/*|*/ 4, 7, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 1, 0,/*|*/ 5, 7, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 1, 0,/*|*/ 6, 7, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 1, 0,/*|*/ 7, 7, 5, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0001'1111}, argument_none{}},
    simd_swizzle_op{{0, 0, 2, 0,/*|*/ 4, 4, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 2, 0,/*|*/ 5, 4, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 2, 0,/*|*/ 6, 4, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 2, 0,/*|*/ 7, 4, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 2, 0,/*|*/ 4, 5, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 2, 0,/*|*/ 5, 5, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 2, 0,/*|*/ 6, 5, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 2, 0,/*|*/ 7, 5, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 2, 0,/*|*/ 4, 6, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 2, 0,/*|*/ 5, 6, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 2, 0,/*|*/ 6, 6, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 2, 0,/*|*/ 7, 6, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 2, 0,/*|*/ 4, 7, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 2, 0,/*|*/ 5, 7, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 2, 0,/*|*/ 6, 7, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 2, 0,/*|*/ 7, 7, 6, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0010'1111}, argument_none{}},
    simd_swizzle_op{{0, 0, 3, 0,/*|*/ 4, 4, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 3, 0,/*|*/ 5, 4, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 3, 0,/*|*/ 6, 4, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 3, 0,/*|*/ 7, 4, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 3, 0,/*|*/ 4, 5, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 3, 0,/*|*/ 5, 5, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 3, 0,/*|*/ 6, 5, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 3, 0,/*|*/ 7, 5, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 3, 0,/*|*/ 4, 6, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 3, 0,/*|*/ 5, 6, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 3, 0,/*|*/ 6, 6, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 3, 0,/*|*/ 7, 6, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 3, 0,/*|*/ 4, 7, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 3, 0,/*|*/ 5, 7, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 3, 0,/*|*/ 6, 7, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 3, 0,/*|*/ 7, 7, 7, 4}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0011'1111}, argument_none{}},

    simd_swizzle_op{{0, 0, 0, 1,/*|*/ 4, 4, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 0, 1,/*|*/ 5, 4, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 0, 1,/*|*/ 6, 4, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 0, 1,/*|*/ 7, 4, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 0, 1,/*|*/ 4, 5, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 0, 1,/*|*/ 5, 5, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 0, 1,/*|*/ 6, 5, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 0, 1,/*|*/ 7, 5, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 0, 1,/*|*/ 4, 6, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 0, 1,/*|*/ 5, 6, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 0, 1,/*|*/ 6, 6, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 0, 1,/*|*/ 7, 6, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 0, 1,/*|*/ 4, 7, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 0, 1,/*|*/ 5, 7, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 0, 1,/*|*/ 6, 7, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 0, 1,/*|*/ 7, 7, 4, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0100'1111}, argument_none{}},
    simd_swizzle_op{{0, 0, 1, 1,/*|*/ 4, 4, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 1, 1,/*|*/ 5, 4, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 1, 1,/*|*/ 6, 4, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 1, 1,/*|*/ 7, 4, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 1, 1,/*|*/ 4, 5, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 1, 1,/*|*/ 5, 5, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 1, 1,/*|*/ 6, 5, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 1, 1,/*|*/ 7, 5, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 1, 1,/*|*/ 4, 6, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 1, 1,/*|*/ 5, 6, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 1, 1,/*|*/ 6, 6, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 1, 1,/*|*/ 7, 6, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 1, 1,/*|*/ 4, 7, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 1, 1,/*|*/ 5, 7, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 1, 1,/*|*/ 6, 7, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 1, 1,/*|*/ 7, 7, 5, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0101'1111}, argument_none{}},
    simd_swizzle_op{{0, 0, 2, 1,/*|*/ 4, 4, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 2, 1,/*|*/ 5, 4, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 2, 1,/*|*/ 6, 4, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 2, 1,/*|*/ 7, 4, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 2, 1,/*|*/ 4, 5, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 2, 1,/*|*/ 5, 5, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 2, 1,/*|*/ 6, 5, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 2, 1,/*|*/ 7, 5, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 2, 1,/*|*/ 4, 6, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 2, 1,/*|*/ 5, 6, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 2, 1,/*|*/ 6, 6, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 2, 1,/*|*/ 7, 6, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 2, 1,/*|*/ 4, 7, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 2, 1,/*|*/ 5, 7, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 2, 1,/*|*/ 6, 7, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 2, 1,/*|*/ 7, 7, 6, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0110'1111}, argument_none{}},
    simd_swizzle_op{{0, 0, 3, 1,/*|*/ 4, 4, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 3, 1,/*|*/ 5, 4, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 3, 1,/*|*/ 6, 4, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 3, 1,/*|*/ 7, 4, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 3, 1,/*|*/ 4, 5, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 3, 1,/*|*/ 5, 5, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 3, 1,/*|*/ 6, 5, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 3, 1,/*|*/ 7, 5, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 3, 1,/*|*/ 4, 6, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 3, 1,/*|*/ 5, 6, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 3, 1,/*|*/ 6, 6, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 3, 1,/*|*/ 7, 6, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 3, 1,/*|*/ 4, 7, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 3, 1,/*|*/ 5, 7, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 3, 1,/*|*/ 6, 7, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 3, 1,/*|*/ 7, 7, 7, 5}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b0111'1111}, argument_none{}},

    simd_swizzle_op{{0, 0, 0, 2,/*|*/ 4, 4, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 0, 2,/*|*/ 5, 4, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 0, 2,/*|*/ 6, 4, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 0, 2,/*|*/ 7, 4, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 0, 2,/*|*/ 4, 5, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 0, 2,/*|*/ 5, 5, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 0, 2,/*|*/ 6, 5, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 0, 2,/*|*/ 7, 5, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 0, 2,/*|*/ 4, 6, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 0, 2,/*|*/ 5, 6, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 0, 2,/*|*/ 6, 6, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 0, 2,/*|*/ 7, 6, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 0, 2,/*|*/ 4, 7, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 0, 2,/*|*/ 5, 7, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 0, 2,/*|*/ 6, 7, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 0, 2,/*|*/ 7, 7, 4, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1000'1111}, argument_none{}},
    simd_swizzle_op{{0, 0, 1, 2,/*|*/ 4, 4, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 1, 2,/*|*/ 5, 4, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 1, 2,/*|*/ 6, 4, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 1, 2,/*|*/ 7, 4, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 1, 2,/*|*/ 4, 5, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 1, 2,/*|*/ 5, 5, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 1, 2,/*|*/ 6, 5, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 1, 2,/*|*/ 7, 5, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 1, 2,/*|*/ 4, 6, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 1, 2,/*|*/ 5, 6, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 1, 2,/*|*/ 6, 6, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 1, 2,/*|*/ 7, 6, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 1, 2,/*|*/ 4, 7, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 1, 2,/*|*/ 5, 7, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 1, 2,/*|*/ 6, 7, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 1, 2,/*|*/ 7, 7, 5, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1001'1111}, argument_none{}},
    simd_swizzle_op{{0, 0, 2, 2,/*|*/ 4, 4, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 2, 2,/*|*/ 5, 4, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 2, 2,/*|*/ 6, 4, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 2, 2,/*|*/ 7, 4, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 2, 2,/*|*/ 4, 5, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 2, 2,/*|*/ 5, 5, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 2, 2,/*|*/ 6, 5, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 2, 2,/*|*/ 7, 5, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 2, 2,/*|*/ 4, 6, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 2, 2,/*|*/ 5, 6, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 2, 2,/*|*/ 6, 6, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 2, 2,/*|*/ 7, 6, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 2, 2,/*|*/ 4, 7, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 2, 2,/*|*/ 5, 7, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 2, 2,/*|*/ 6, 7, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 2, 2,/*|*/ 7, 7, 6, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1010'1111}, argument_none{}},
    simd_swizzle_op{{0, 0, 3, 2,/*|*/ 4, 4, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 3, 2,/*|*/ 5, 4, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 3, 2,/*|*/ 6, 4, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 3, 2,/*|*/ 7, 4, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 3, 2,/*|*/ 4, 5, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 3, 2,/*|*/ 5, 5, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 3, 2,/*|*/ 6, 5, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 3, 2,/*|*/ 7, 5, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 3, 2,/*|*/ 4, 6, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 3, 2,/*|*/ 5, 6, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 3, 2,/*|*/ 6, 6, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 3, 2,/*|*/ 7, 6, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 3, 2,/*|*/ 4, 7, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 3, 2,/*|*/ 5, 7, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 3, 2,/*|*/ 6, 7, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 3, 2,/*|*/ 7, 7, 7, 6}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1011'1111}, argument_none{}},

    simd_swizzle_op{{0, 0, 0, 3,/*|*/ 4, 4, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 0, 3,/*|*/ 5, 4, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 0, 3,/*|*/ 6, 4, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 0, 3,/*|*/ 7, 4, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 0, 3,/*|*/ 4, 5, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 0, 3,/*|*/ 5, 5, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 0, 3,/*|*/ 6, 5, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 0, 3,/*|*/ 7, 5, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 0, 3,/*|*/ 4, 6, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 0, 3,/*|*/ 5, 6, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 0, 3,/*|*/ 6, 6, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 0, 3,/*|*/ 7, 6, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 0, 3,/*|*/ 4, 7, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 0, 3,/*|*/ 5, 7, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 0, 3,/*|*/ 6, 7, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 0, 3,/*|*/ 7, 7, 4, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1100'1111}, argument_none{}},
    simd_swizzle_op{{0, 0, 1, 3,/*|*/ 4, 4, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 1, 3,/*|*/ 5, 4, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 1, 3,/*|*/ 6, 4, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 1, 3,/*|*/ 7, 4, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 1, 3,/*|*/ 4, 5, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 1, 3,/*|*/ 5, 5, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 1, 3,/*|*/ 6, 5, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 1, 3,/*|*/ 7, 5, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 1, 3,/*|*/ 4, 6, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 1, 3,/*|*/ 5, 6, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 1, 3,/*|*/ 6, 6, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 1, 3,/*|*/ 7, 6, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 1, 3,/*|*/ 4, 7, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 1, 3,/*|*/ 5, 7, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 1, 3,/*|*/ 6, 7, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 1, 3,/*|*/ 7, 7, 5, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1101'1111}, argument_none{}},
    simd_swizzle_op{{0, 0, 2, 3,/*|*/ 4, 4, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 2, 3,/*|*/ 5, 4, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 2, 3,/*|*/ 6, 4, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 2, 3,/*|*/ 7, 4, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 2, 3,/*|*/ 4, 5, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 2, 3,/*|*/ 5, 5, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 2, 3,/*|*/ 6, 5, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 2, 3,/*|*/ 7, 5, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 2, 3,/*|*/ 4, 6, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 2, 3,/*|*/ 5, 6, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 2, 3,/*|*/ 6, 6, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 2, 3,/*|*/ 7, 6, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 2, 3,/*|*/ 4, 7, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 2, 3,/*|*/ 5, 7, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 2, 3,/*|*/ 6, 7, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 2, 3,/*|*/ 7, 7, 6, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1110'1111}, argument_none{}},
    simd_swizzle_op{{0, 0, 3, 3,/*|*/ 4, 4, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'0000}, argument_none{}},
    simd_swizzle_op{{1, 0, 3, 3,/*|*/ 5, 4, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'0001}, argument_none{}},
    simd_swizzle_op{{2, 0, 3, 3,/*|*/ 6, 4, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'0010}, argument_none{}},
    simd_swizzle_op{{3, 0, 3, 3,/*|*/ 7, 4, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'0011}, argument_none{}},
    simd_swizzle_op{{0, 1, 3, 3,/*|*/ 4, 5, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'0100}, argument_none{}},
    simd_swizzle_op{{1, 1, 3, 3,/*|*/ 5, 5, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'0101}, argument_none{}},
    simd_swizzle_op{{2, 1, 3, 3,/*|*/ 6, 5, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'0110}, argument_none{}},
    simd_swizzle_op{{3, 1, 3, 3,/*|*/ 7, 5, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'0111}, argument_none{}},
    simd_swizzle_op{{0, 2, 3, 3,/*|*/ 4, 6, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'1000}, argument_none{}},
    simd_swizzle_op{{1, 2, 3, 3,/*|*/ 5, 6, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'1001}, argument_none{}},
    simd_swizzle_op{{2, 2, 3, 3,/*|*/ 6, 6, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'1010}, argument_none{}},
    simd_swizzle_op{{3, 2, 3, 3,/*|*/ 7, 6, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'1011}, argument_none{}},
    simd_swizzle_op{{0, 3, 3, 3,/*|*/ 4, 7, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'1100}, argument_none{}},
    simd_swizzle_op{{1, 3, 3, 3,/*|*/ 5, 7, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'1101}, argument_none{}},
    simd_swizzle_op{{2, 3, 3, 3,/*|*/ 6, 7, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'1110}, argument_none{}},
    simd_swizzle_op{{3, 3, 3, 3,/*|*/ 7, 7, 7, 7}, "_mm256_shuffle_epi32", argument_m256i{}, argument_imm8{0b1111'1111}, argument_none{}},
};

template <size_t> struct simd_vars;
template <size_t count> struct simd_swizzle_solver;

template <template <size_t> typename simd_swizzle_solver_t = simd_swizzle_solver, template <size_t> typename simd_vars_t = simd_vars>
void simd_swizzle_solver_shuffle_epi32_test(auto const & transitions)
{
    // _mm256_shuffle_epi32_0b0000'0000
    {
        auto && [memory, final_state] = simd_swizzle_solver_t<1>{}.find(
            transitions,
            simd_vars_t<1>::iota(), // start state
            simd_vars_t<1>{0, 0, 0, 0,/*|*/ 4, 4, 4, 4}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_imm8, argument_none>>(final_state.ops[0]);
        assert(state0.name == "_mm256_shuffle_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg1.constant == 0b0000'0000);
        assert((memory == simd_vars_t<1>{0, 0, 0, 0,/*|*/ 4, 4, 4, 4}));
    }

    // _mm256_shuffle_epi32_0b0111'0101
    {
        auto && [memory, final_state] = simd_swizzle_solver_t<1>{}.find(
            transitions,
            simd_vars_t<1>::iota(), // start state
            simd_vars_t<1>{1, 1, 3, 1,/*|*/ 5, 5, 7, 5}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_imm8, argument_none>>(final_state.ops[0]);
        assert(state0.name == "_mm256_shuffle_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg1.constant == 0b0111'0101);
        assert((memory == simd_vars_t<1>{1, 1, 3, 1,/*|*/ 5, 5, 7, 5}));
    }

    // _mm256_shuffle_epi32_0b0100'1111
    {
        auto && [memory, final_state] = simd_swizzle_solver_t<1>{}.find(
            transitions,
            simd_vars_t<1>::iota(), // start state
            simd_vars_t<1>{3, 3, 0, 1,/*|*/ 7, 7, 4, 5}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_imm8, argument_none>>(final_state.ops[0]);
        assert(state0.name == "_mm256_shuffle_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg1.constant == 0b0100'1111);
        assert((memory == simd_vars_t<1>{3, 3, 0, 1,/*|*/ 7, 7, 4, 5}));
    }
}
