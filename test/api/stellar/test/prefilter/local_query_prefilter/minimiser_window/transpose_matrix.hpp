
#pragma once

#include <span>

#include "diagnostics.hpp"
#include "simd_fn.hpp"

void transpose_matrix_32x4x4_omp(std::span<int32x4_t, 4> matrix)
{
    #pragma omp simd
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            // std::cout << "M[" << i << "," << j << "]: " << matrix[i][j];
            // std::cout << " <-> M[" << j << "," << i << "]: " << matrix[j][i] << std::endl;
            int32_t tmp1 = matrix[i][j];
            int32_t tmp2 = matrix[j][i];
            matrix[i][j] = tmp2;
            matrix[j][i] = tmp1;
        }
    }
}

void transpose_matrix_32x4x4_sse4(std::span<int32x4_t, 4> matrix)
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

void transpose_matrix_32x4x4_avx2(std::span<int32x4_t, 4> matrix)
{
    // std::cout << "matrix: " << std::endl;
    // print_simd(matrix[0], 4);
    // print_simd(matrix[1], 4);
    // print_simd(matrix[2], 4);
    // print_simd(matrix[3], 4);
    __m256i * matrix_data = reinterpret_cast<__m256i *>(matrix.data());
    __m256i row01 = _mm256_loadu_si256(matrix_data);
    __m256i row23 = _mm256_loadu_si256(matrix_data + 1);
    // __m256i row01 = _mm256_set_m128i((__m128i)matrix[1], (__m128i)matrix[0]);
    // __m256i row23 = _mm256_set_m128i((__m128i)matrix[3], (__m128i)matrix[2]);
    // print_simd(row01, 8);
    // print_simd(row23, 8);

    __m256i tmp01 = _mm256_unpacklo_epi32(row01, row23);
    __m256i tmp23 = _mm256_unpackhi_epi32(row01, row23);
    // std::cout << "tmp01/tmp23: " << std::endl;
    // print_simd(tmp01, 8);
    // print_simd(tmp23, 8);

    __m256i shuffleCells = (__m256i)int32x8_t{0, 4, 1, 5, 2, 6, 3, 7};
    row01 = _mm256_permutevar8x32_epi32(tmp01, shuffleCells);
    row23 = _mm256_permutevar8x32_epi32(tmp23, shuffleCells);

    // row01 = _mm256_shuffle_epi32(tmp01, 0b0000'0011);
    // std::cout << "row01/row23: " << std::endl;
    // print_simd(row01, 8);
    // print_simd(row23, 8);
    // // __m128 tmp3, tmp2, tmp1, tmp0;
    // tmp0 = _mm_unpacklo_ps((__m128)matrix[0], (__m128)matrix[1]);
    // tmp2 = _mm_unpacklo_ps((__m128)matrix[2], (__m128)matrix[3]);
    // tmp1 = _mm_unpackhi_ps((__m128)matrix[0], (__m128)matrix[1]);
    // tmp3 = _mm_unpackhi_ps((__m128)matrix[2], (__m128)matrix[3]);
    // matrix[0] = (int32x4_t)_mm256_extractf128_si256(row01, 0);
    // matrix[1] = (int32x4_t)_mm256_extractf128_si256(row01, 1);
    // matrix[2] = (int32x4_t)_mm256_extractf128_si256(row23, 0);
    // matrix[3] = (int32x4_t)_mm256_extractf128_si256(row23, 1);
    _mm256_storeu_si256(matrix_data, row01);
    _mm256_storeu_si256(matrix_data + 1, row23);
    // std::cout << "matrix: " << std::endl;
    // print_simd(matrix[0], 4);
    // print_simd(matrix[1], 4);
    // print_simd(matrix[2], 4);
    // print_simd(matrix[3], 4);
}

void transpose_matrix_32x4x4_avx2_gather(std::span<int32x4_t, 4> matrix)
{
    int32_t * matrix_data = reinterpret_cast<int32_t *>(matrix.data());

    __m256i new_row0 = _mm256_i32gather_epi32(matrix_data, (__m256i)(int32x8_t{0, 4, 8, 12, 1, 5, 9, 13}), sizeof(int32_t));
    __m256i new_row1 = _mm256_i32gather_epi32(matrix_data, (__m256i)(int32x8_t{2, 6, 10, 14, 3, 7, 11, 15}), sizeof(int32_t));

    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data), new_row0);
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data) + 1, new_row1);
}

void transpose_matrix_32x8x8_avx2(std::span<int32x8_t, 8> matrix)
{
    // std::cout << "matrix: " << std::endl;
    // print_simd(matrix[0], 4);
    // print_simd(matrix[1], 4);
    // print_simd(matrix[2], 4);
    // print_simd(matrix[3], 4);
    __m256i * matrix_data = reinterpret_cast<__m256i *>(matrix.data());
    __m256i row0 = _mm256_loadu_si256(matrix_data);
    __m256i row1 = _mm256_loadu_si256(matrix_data + 1);
    __m256i row2 = _mm256_loadu_si256(matrix_data + 2);
    __m256i row3 = _mm256_loadu_si256(matrix_data + 3);
    __m256i row4 = _mm256_loadu_si256(matrix_data + 4);
    __m256i row5 = _mm256_loadu_si256(matrix_data + 5);
    __m256i row6 = _mm256_loadu_si256(matrix_data + 6);
    __m256i row7 = _mm256_loadu_si256(matrix_data + 7);

    __m256i a0 = _mm256_permute2f128_si256(row0, row4, 0b0010'0000);
    __m256i a1 = _mm256_permute2f128_si256(row1, row5, 0b0010'0000);
    __m256i a2 = _mm256_permute2f128_si256(row2, row6, 0b0010'0000);
    __m256i a3 = _mm256_permute2f128_si256(row3, row7, 0b0010'0000);

    __m256i b02lo = _mm256_unpacklo_epi32(a0, a2);
    __m256i b13lo = _mm256_unpacklo_epi32(a1, a3);
    __m256i new_row0 = _mm256_unpacklo_epi32(b02lo, b13lo);
    __m256i new_row1 = _mm256_unpackhi_epi32(b02lo, b13lo);

    __m256i b02hi = _mm256_unpackhi_epi32(a0, a2);
    __m256i b13hi = _mm256_unpackhi_epi32(a1, a3);
    __m256i new_row2 = _mm256_unpacklo_epi32(b02hi, b13hi);
    __m256i new_row3 = _mm256_unpackhi_epi32(b02hi, b13hi);

    __m256i a4 = _mm256_permute2f128_si256(row0, row4, 0b0011'0001);
    __m256i a5 = _mm256_permute2f128_si256(row1, row5, 0b0011'0001);
    __m256i a6 = _mm256_permute2f128_si256(row2, row6, 0b0011'0001);
    __m256i a7 = _mm256_permute2f128_si256(row3, row7, 0b0011'0001);

    __m256i b46lo = _mm256_unpacklo_epi32(a4, a6);
    __m256i b57lo = _mm256_unpacklo_epi32(a5, a7);
    __m256i new_row4 = _mm256_unpacklo_epi32(b46lo, b57lo);
    __m256i new_row5 = _mm256_unpackhi_epi32(b46lo, b57lo);

    __m256i b46hi = _mm256_unpackhi_epi32(a4, a6);
    __m256i b57hi = _mm256_unpackhi_epi32(a5, a7);
    __m256i new_row6 = _mm256_unpacklo_epi32(b46hi, b57hi);
    __m256i new_row7 = _mm256_unpackhi_epi32(b46hi, b57hi);

    _mm256_storeu_si256(matrix_data, new_row0);
    _mm256_storeu_si256(matrix_data + 1, new_row1);
    _mm256_storeu_si256(matrix_data + 2, new_row2);
    _mm256_storeu_si256(matrix_data + 3, new_row3);
    _mm256_storeu_si256(matrix_data + 4, new_row4);
    _mm256_storeu_si256(matrix_data + 5, new_row5);
    _mm256_storeu_si256(matrix_data + 6, new_row6);
    _mm256_storeu_si256(matrix_data + 7, new_row7);

    // print_simd(b02hi, 8);
    // print_simd(b13hi, 8);
    // print_simd(new_row0, 8);
    // print_simd(new_row1, 8);
    // print_simd(new_row2, 8);
    // print_simd(new_row3, 8);
    // print_simd(new_row4, 8);
    // print_simd(new_row5, 8);
    // __m256i row01 = _mm256_set_m128i((__m128i)matrix[1], (__m128i)matrix[0]);
    // __m256i row23 = _mm256_set_m128i((__m128i)matrix[3], (__m128i)matrix[2]);
    // print_simd(row23, 8);

    // __m256i tmp01 = _mm256_unpacklo_epi32(row01, row23);
    // __m256i tmp23 = _mm256_unpackhi_epi32(row01, row23);
    // // std::cout << "tmp01/tmp23: " << std::endl;
    // // print_simd(tmp01, 8);
    // // print_simd(tmp23, 8);
    //
    // __m256i shuffleCells = (__m256i)int32x8_t{0, 4, 1, 5, 2, 6, 3, 7};
    // row01 = _mm256_permutevar8x32_epi32(tmp01, shuffleCells);
    // row23 = _mm256_permutevar8x32_epi32(tmp23, shuffleCells);

    // row01 = _mm256_shuffle_epi32(tmp01, 0b0000'0011);
    // std::cout << "row01/row23: " << std::endl;
    // print_simd(row01, 8);
    // print_simd(row23, 8);
    // // __m128 tmp3, tmp2, tmp1, tmp0;
    // tmp0 = _mm_unpacklo_ps((__m128)matrix[0], (__m128)matrix[1]);
    // tmp2 = _mm_unpacklo_ps((__m128)matrix[2], (__m128)matrix[3]);
    // tmp1 = _mm_unpackhi_ps((__m128)matrix[0], (__m128)matrix[1]);
    // tmp3 = _mm_unpackhi_ps((__m128)matrix[2], (__m128)matrix[3]);
    // matrix[0] = (int32x4_t)_mm256_extractf128_si256(row01, 0);
    // matrix[1] = (int32x4_t)_mm256_extractf128_si256(row01, 1);
    // matrix[2] = (int32x4_t)_mm256_extractf128_si256(row23, 0);
    // matrix[3] = (int32x4_t)_mm256_extractf128_si256(row23, 1);
    // _mm256_storeu_si256(matrix_data, row01);
    // _mm256_storeu_si256(matrix_data + 1, row23);
    // std::cout << "matrix: " << std::endl;
    // print_simd(matrix[0], 4);
    // print_simd(matrix[1], 4);
    // print_simd(matrix[2], 4);
    // print_simd(matrix[3], 4);
}

void transpose_matrix_32x8x8_avx2_gather(std::span<int32x8_t, 8> matrix)
{
    int32_t * matrix_data = reinterpret_cast<int32_t *>(matrix.data());

    int32x8_t offset = int32x8_t{0, 1, 2, 3, 4, 5, 6, 7} * 8;

    __m256i new_row0 = _mm256_i32gather_epi32(matrix_data, (__m256i)(offset + 0), sizeof(int32_t));
    __m256i new_row1 = _mm256_i32gather_epi32(matrix_data, (__m256i)(offset + 1), sizeof(int32_t));
    __m256i new_row2 = _mm256_i32gather_epi32(matrix_data, (__m256i)(offset + 2), sizeof(int32_t));
    __m256i new_row3 = _mm256_i32gather_epi32(matrix_data, (__m256i)(offset + 3), sizeof(int32_t));
    __m256i new_row4 = _mm256_i32gather_epi32(matrix_data, (__m256i)(offset + 4), sizeof(int32_t));
    __m256i new_row5 = _mm256_i32gather_epi32(matrix_data, (__m256i)(offset + 5), sizeof(int32_t));
    __m256i new_row6 = _mm256_i32gather_epi32(matrix_data, (__m256i)(offset + 6), sizeof(int32_t));
    __m256i new_row7 = _mm256_i32gather_epi32(matrix_data, (__m256i)(offset + 7), sizeof(int32_t));

    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data) + 0, new_row0);
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data) + 1, new_row1);
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data) + 2, new_row2);
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data) + 3, new_row3);
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data) + 4, new_row4);
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data) + 5, new_row5);
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data) + 6, new_row6);
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data) + 7, new_row7);
}

void transpose_matrix_32x8x8_omp(std::span<int32x8_t, 8> matrix)
{
    #pragma omp simd
    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            int32_t tmp1 = matrix[i][j];
            int32_t tmp2 = matrix[j][i];
            matrix[i][j] = tmp2;
            matrix[j][i] = tmp1;
        }
    }
}

void transpose_matrix_32x6x4_omp(std::span<int32x4_t, 6> matrix)
{
    // for (int i = 0; i < 6; ++i)
    // {
    //     for (int j = 0; j < 4; ++j)
    //     {
    //         std::cout << "M[" << i << "," << j << "]: " << matrix[i][j] << std::endl;
    //     }
    //     std::cout << std::endl;
    // }
    int * matrix_ptr = reinterpret_cast<int *>(matrix.data());
    std::array<int, 6 * 4> matrix_cpy;
    std::copy(matrix_ptr, matrix_ptr + 6 * 4, matrix_cpy.data());

    // for (int i = 0; i < matrix_cpy.size(); ++i)
    // {
    //     std::cout << "Mc[" << i << "]: " << matrix_cpy[i] << std::endl;
    // }
    // std::cout << std::endl;

    #pragma omp simd
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            matrix[i][j] = matrix_cpy[j * 6 + i];
            // std::cout << "M[" << i << "," << j << "]: " << matrix_cpy[j * 6 + i] << std::endl;
        }
    }
}

void transpose_matrix_32x6x8_omp(std::span<int32x8_t, 6> matrix)
{
    int * matrix_ptr = reinterpret_cast<int *>(matrix.data());
    std::array<int, 6 * 8> matrix_cpy;
    std::copy(matrix_ptr, matrix_ptr + 6 * 8, matrix_cpy.data());

    #pragma omp simd
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 8; ++j)
        {
            matrix[i][j] = matrix_cpy[j * 6 + i];
        }
    }
}

void transpose_matrix_32x6x8_avx2(std::span<int32x8_t, 6> matrix)
{
    __m256i * matrix_data = reinterpret_cast<__m256i *>(matrix.data());
    __m256i row0 = _mm256_loadu_si256(matrix_data);
    __m256i row1 = _mm256_loadu_si256(matrix_data + 1);
    __m256i row2 = _mm256_loadu_si256(matrix_data + 2);
    __m256i row3 = _mm256_loadu_si256(matrix_data + 3);
    __m256i row4 = _mm256_loadu_si256(matrix_data + 4);
    __m256i row5 = _mm256_loadu_si256(matrix_data + 5);

    //  0,  1,  2,  3,  4,  5,
    //  6,  7,  8,  9, 10, 11,
    // 12, 13, 14, 15, 16, 17,
    // 18, 19, 20, 21, 22, 23,
    // 24, 25, 26, 27, 28, 29,
    // 30, 31, 32, 33, 34, 35,
    // 36, 37, 38, 39, 40, 41,
    // 42, 43, 44, 45, 46, 47,

    //* 0,  1,  2,  3,  4,  5,* 6,  7,
    //  8,  9, 10, 11,*12, 13, 14, 15,
    // 16, 17,*18, 19, 20, 21, 22, 23,
    //*24, 25, 26, 27, 28, 29,*30, 31,
    // 32, 33, 34, 35,*36, 37, 38, 39,
    // 40, 41,*42, 43, 44, 45, 46, 47,

    /*0:*/ __m256i s03_lo = _mm256_unpacklo_epi32(row0, row3); //* 0,*24,  1, 25,  4, 28,  5, 29,
    /*1:*/ __m256i s03_hi = _mm256_unpackhi_epi32(row0, row3); //  2, 26,  3, 27,* 6,*30,  7, 31,
    /*2:*/ __m256i s14_lo = _mm256_unpacklo_epi32(row1, row4); //  8, 32,  9, 33,*12,*36, 13, 37,
    /*3:*/ __m256i s14_hi = _mm256_unpackhi_epi32(row1, row4); // 10, 34, 11, 35, 14, 38, 15, 39,
    /*4:*/ __m256i s25_lo = _mm256_unpacklo_epi32(row2, row5); // 16, 40, 17, 41, 20, 44, 21, 45,
    /*5:*/ __m256i s25_hi = _mm256_unpackhi_epi32(row2, row5); //*18,*42, 19, 43, 22, 46, 23, 47,

    __m256i t05_lo = _mm256_unpacklo_epi32(s03_lo, s25_hi); // *0,*18,*24,*42,  4, 22, 28, 46,
    __m256i t05_hi = _mm256_unpackhi_epi32(s03_lo, s25_hi); //  1, 19, 25, 43,  5, 23, 29, 47,
    __m256i t12_lo = _mm256_unpacklo_epi32(s03_hi, s14_lo); //  2,  8, 26, 32, *6,*12,*30,*36,
    __m256i t12_hi = _mm256_unpackhi_epi32(s03_hi, s14_lo); //  3,  9, 27, 33,  7, 13, 31, 37,
    __m256i t34_lo = _mm256_unpacklo_epi32(s14_hi, s25_lo); // 10, 16, 34, 40, 14, 20, 38, 44,
    __m256i t34_hi = _mm256_unpackhi_epi32(s14_hi, s25_lo); // 11, 17, 35, 41, 15, 21, 39, 45,

    __m256i unsorted_row0 = _mm256_blend_epi32(t05_lo, t12_lo, 0b11110000); // *0,*18,*24,*42, *6,*12,*30,*36,
    __m256i unsorted_row1 = _mm256_blend_epi32(t05_hi, t12_hi, 0b11110000); //  1, 19, 25, 43,  7, 13, 31, 37
    __m256i unsorted_row2 = _mm256_blend_epi32(t12_lo, t34_lo, 0b11110000); //  2,  8, 26, 32, 14, 20, 38, 44,
    __m256i unsorted_row3 = _mm256_blend_epi32(t12_hi, t34_hi, 0b11110000); //  3,  9, 27, 33, 15, 21, 39, 45,
    __m256i unsorted_row4 = _mm256_blend_epi32(t34_lo, t05_lo, 0b11110000); // 10, 16, 34, 40,  4, 22, 28, 46,
    __m256i unsorted_row5 = _mm256_blend_epi32(t34_hi, t05_hi, 0b11110000); // 11, 17, 35, 41,  5, 23, 29, 47,

    __m256i new_row0 = _mm256_permutevar8x32_epi32(unsorted_row0, (__m256i)(int32x8_t{0, 4, 5, 1, 2, 6, 7, 3}));
    __m256i new_row1 = _mm256_permutevar8x32_epi32(unsorted_row1, (__m256i)(int32x8_t{0, 4, 5, 1, 2, 6, 7, 3}));
    __m256i new_row2 = _mm256_permutevar8x32_epi32(unsorted_row2, (__m256i)(int32x8_t{0, 1, 4, 5, 2, 3, 6, 7}));
    __m256i new_row3 = _mm256_permutevar8x32_epi32(unsorted_row3, (__m256i)(int32x8_t{0, 1, 4, 5, 2, 3, 6, 7}));
    __m256i new_row4 = _mm256_permutevar8x32_epi32(unsorted_row4, (__m256i)(int32x8_t{4, 0, 1, 5, 6, 2, 3, 7}));
    __m256i new_row5 = _mm256_permutevar8x32_epi32(unsorted_row5, (__m256i)(int32x8_t{4, 0, 1, 5, 6, 2, 3, 7}));

    _mm256_storeu_si256(matrix_data + 0, new_row0);
    _mm256_storeu_si256(matrix_data + 1, new_row1);
    _mm256_storeu_si256(matrix_data + 2, new_row2);
    _mm256_storeu_si256(matrix_data + 3, new_row3);
    _mm256_storeu_si256(matrix_data + 4, new_row4);
    _mm256_storeu_si256(matrix_data + 5, new_row5);
}

void transpose_matrix_32x6x4_sse4(std::span<int32x4_t, 6> matrix)
{
    __m128 s03_lo = _mm_unpacklo_ps((__m128)matrix[0], (__m128)matrix[3]);
    __m128 s14_hi = _mm_unpackhi_ps((__m128)matrix[1], (__m128)matrix[4]);

    __m128 s03_hi = _mm_unpackhi_ps((__m128)matrix[0], (__m128)matrix[3]);
    __m128 s25_lo = _mm_unpacklo_ps((__m128)matrix[2], (__m128)matrix[5]);

    __m128 s14_lo = _mm_unpacklo_ps((__m128)matrix[1], (__m128)matrix[4]);
    __m128 s25_hi = _mm_unpackhi_ps((__m128)matrix[2], (__m128)matrix[5]);

    matrix[0] = (int32x4_t)_mm_unpacklo_ps(s03_lo, s14_hi);
    matrix[1] = (int32x4_t)_mm_unpackhi_ps(s03_lo, s14_hi);
    matrix[2] = (int32x4_t)_mm_unpacklo_ps(s03_hi, s25_lo);
    matrix[3] = (int32x4_t)_mm_unpackhi_ps(s03_hi, s25_lo);
    matrix[4] = (int32x4_t)_mm_unpacklo_ps(s14_lo, s25_hi);
    matrix[5] = (int32x4_t)_mm_unpackhi_ps(s14_lo, s25_hi);
}

void transpose_matrix_32x6x4_avx2(std::span<int32x4_t, 6> matrix)
{
    __m256i * matrix_data = reinterpret_cast<__m256i *>(matrix.data());
    __m256i row0 = _mm256_loadu_si256(matrix_data);
    __m256i row1 = _mm256_loadu_si256(matrix_data + 1);
    __m256i row2 = _mm256_loadu_si256(matrix_data + 2);

    __m256i s0 = _mm256_blend_epi32(row0, row1, 0b0000'1100);
    __m256i s1 = _mm256_blend_epi32(row1, row2, 0b0000'1100);
    __m256i s2 = _mm256_blend_epi32(row2, row0, 0b0000'1100);

    __m256i t0 = _mm256_blend_epi32(s0, s1, 0b0011'1100);
    __m256i t1 = _mm256_blend_epi32(s1, s2, 0b0011'1100);
    __m256i t2 = _mm256_blend_epi32(s2, s0, 0b0011'1100);

    __m256i new_row0 = _mm256_permutevar8x32_epi32(t0, (__m256i)int32x8_t{0, 6, 4, 2, 1, 7, 5, 3});
    __m256i new_row1 = _mm256_permutevar8x32_epi32(t1, (__m256i)int32x8_t{2, 0, 6, 4, 3, 1, 7, 5});
    __m256i new_row2 = _mm256_permutevar8x32_epi32(t2, (__m256i)int32x8_t{4, 2, 0, 6, 5, 3, 1, 7});

    _mm256_storeu_si256(matrix_data, new_row0);
    _mm256_storeu_si256(matrix_data + 1, new_row1);
    _mm256_storeu_si256(matrix_data + 2, new_row2);
}

void transpose_matrix_32x6x4_avx2_gather(std::span<int32x4_t, 6> matrix)
{
    int32_t * matrix_data = reinterpret_cast<int32_t *>(matrix.data());

    __m256i new_row0 = _mm256_i32gather_epi32(matrix_data, (__m256i)(int32x8_t{0, 6, 12, 18, 1, 7, 13, 19}), sizeof(int32_t));
    __m256i new_row1 = _mm256_i32gather_epi32(matrix_data, (__m256i)(int32x8_t{2, 8, 14, 20, 3, 9, 15, 21}), sizeof(int32_t));
    __m256i new_row2 = _mm256_i32gather_epi32(matrix_data, (__m256i)(int32x8_t{4,10, 16, 22, 5,11, 17, 23}), sizeof(int32_t));

    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data) + 0, new_row0);
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data) + 1, new_row1);
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(matrix_data) + 2, new_row2);
}

template <auto transpose_matrix_fn>
void transpose_matrix_int32x6x4_test()
{
    {
        std::array<int32x4_t, 6> matrix{};
        transpose_matrix_fn(matrix);
        assert(simd_equal(matrix[0], int32x4_t{}));
        assert(simd_equal(matrix[1], int32x4_t{}));
        assert(simd_equal(matrix[2], int32x4_t{}));
        assert(simd_equal(matrix[3], int32x4_t{}));
        assert(simd_equal(matrix[4], int32x4_t{}));
        assert(simd_equal(matrix[5], int32x4_t{}));
    }

    {
        std::array<int, 6 * 4> matrix
        {
             0,  1,  2,  3,  4,  5,
             6,  7,  8,  9, 10, 11,
            12, 13, 14, 15, 16, 17,
            18, 19, 20, 21, 22, 23
        };

        int32x4_t * matrix_ptr = reinterpret_cast<int32x4_t *>(matrix.data());
        std::span<int32x4_t, 6> matrix_span{matrix_ptr, 6};

        transpose_matrix_fn(matrix_span);
        assert(simd_equal(matrix_span[0], int32x4_t{0, 6, 12, 18}));
        assert(simd_equal(matrix_span[1], int32x4_t{1, 7, 13, 19}));
        assert(simd_equal(matrix_span[2], int32x4_t{2, 8, 14, 20}));
        assert(simd_equal(matrix_span[3], int32x4_t{3, 9, 15, 21}));
        assert(simd_equal(matrix_span[4], int32x4_t{4,10, 16, 22}));
        assert(simd_equal(matrix_span[5], int32x4_t{5,11, 17, 23}));
    }
}

template <auto transpose_matrix_fn>
void transpose_matrix_int32x4x4_test()
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

template <auto transpose_matrix_fn>
void transpose_matrix_int32x6x8_test()
{
    {
        std::array<int32x8_t, 6> matrix{};
        transpose_matrix_fn(matrix);
        assert(simd_equal(matrix[0], int32x8_t{}));
        assert(simd_equal(matrix[1], int32x8_t{}));
        assert(simd_equal(matrix[2], int32x8_t{}));
        assert(simd_equal(matrix[3], int32x8_t{}));
        assert(simd_equal(matrix[4], int32x8_t{}));
        assert(simd_equal(matrix[5], int32x8_t{}));
    }

    {
        std::array<int, 6 * 8> matrix
        {
             0,  1,  2,  3,  4,  5,
             6,  7,  8,  9, 10, 11,
            12, 13, 14, 15, 16, 17,
            18, 19, 20, 21, 22, 23,

            24, 25, 26, 27, 28, 29,
            30, 31, 32, 33, 34, 35,
            36, 37, 38, 39, 40, 41,
            42, 43, 44, 45, 46, 47
        };

        int32x8_t * matrix_ptr = reinterpret_cast<int32x8_t *>(matrix.data());
        std::span<int32x8_t, 6> matrix_span{matrix_ptr, 6};

        transpose_matrix_fn(matrix_span);
        assert(simd_equal(matrix_span[0], int32x8_t{0, 6, 12, 18, 24, 30, 36, 42}));
        assert(simd_equal(matrix_span[1], int32x8_t{1, 7, 13, 19, 25, 31, 37, 43}));
        assert(simd_equal(matrix_span[2], int32x8_t{2, 8, 14, 20, 26, 32, 38, 44}));
        assert(simd_equal(matrix_span[3], int32x8_t{3, 9, 15, 21, 27, 33, 39, 45}));
        assert(simd_equal(matrix_span[4], int32x8_t{4,10, 16, 22, 28, 34, 40, 46}));
        assert(simd_equal(matrix_span[5], int32x8_t{5,11, 17, 23, 29, 35, 41, 47}));
    }
}

template <auto transpose_matrix_fn>
void transpose_matrix_int32x8x8_test()
{
    {
        std::array<int32x8_t, 8> matrix{};
        transpose_matrix_fn(matrix);
        assert(simd_equal(matrix[0], int32x8_t{}));
        assert(simd_equal(matrix[1], int32x8_t{}));
        assert(simd_equal(matrix[2], int32x8_t{}));
        assert(simd_equal(matrix[3], int32x8_t{}));
        assert(simd_equal(matrix[4], int32x8_t{}));
        assert(simd_equal(matrix[5], int32x8_t{}));
        assert(simd_equal(matrix[6], int32x8_t{}));
        assert(simd_equal(matrix[7], int32x8_t{}));
    }

    {
        std::array<int32x8_t, 8> matrix
        {
            int32x8_t{ 0,  1,  2,  3,  4,  5,  6,  7},
            int32x8_t{ 8,  9, 10, 11, 12, 13, 14, 15},
            int32x8_t{16, 17, 18, 19, 20, 21, 22, 23},
            int32x8_t{24, 25, 26, 27, 28, 29, 30, 31},
            int32x8_t{32, 33, 34, 35, 36, 37, 38, 39},
            int32x8_t{40, 41, 42, 43, 44, 45, 46, 47},
            int32x8_t{48, 49, 50, 51, 52, 53, 54, 55},
            int32x8_t{56, 57, 58, 59, 60, 61, 62, 63}
        };
        transpose_matrix_fn(matrix);
        assert(simd_equal(matrix[0], int32x8_t{0, 8, 16, 24, 32, 40, 48, 56}));
        assert(simd_equal(matrix[1], int32x8_t{1, 9, 17, 25, 33, 41, 49, 57}));
        assert(simd_equal(matrix[2], int32x8_t{2, 10, 18, 26, 34, 42, 50, 58}));
        assert(simd_equal(matrix[3], int32x8_t{3, 11, 19, 27, 35, 43, 51, 59}));
        assert(simd_equal(matrix[4], int32x8_t{4, 12, 20, 28, 36, 44, 52, 60}));
        assert(simd_equal(matrix[5], int32x8_t{5, 13, 21, 29, 37, 45, 53, 61}));
        assert(simd_equal(matrix[6], int32x8_t{6, 14, 22, 30, 38, 46, 54, 62}));
        assert(simd_equal(matrix[7], int32x8_t{7, 15, 23, 31, 39, 47, 55, 63}));
    }
}
