
#pragma once
#include <vector> // TODO: remove me, only used for gather_test

#include <cassert>
#include <cmath>
#include <limits>

#include "simd_fn.hpp"

namespace stellar::test
{

inline void forward_kernel_simd(int32x8_t const source, int32x8_t offset, int32x8_t & current_minimiser, int32x8_t & current_offset)
{
    int32x8_mask_t smaller = source < current_minimiser;
    current_minimiser = smaller ? source : current_minimiser;
    current_offset = smaller ? offset : current_offset;
}

void compute_forward_simd(size_t const window_size, int const * __restrict source_ptr, int * __restrict target_ptr, size_t const target_size, int * __restrict offset_ptr)
{
    int32x8_t current_minimiser = simd_set<int32x8_t>(std::numeric_limits<int>::max());
    int32x8_t current_offset = simd_set<int32x8_t>(window_size);

    int32x8_t temporary_minimiser_data[8];
    int32x8_t temporary_offset_data[8];

    for (size_t i = window_size - 1; i > 0; )
    {
        --i;

        // std::cout << "~~~~~~ " << i << std::endl;
        int32x8_t source_column_i = simd_gather(window_size, i, source_ptr);
        // std::cout << "source_column_i: "; print_simd(source_column_i);
        int32x8_t offset = simd_set<int32x8_t>(i);
        // std::cout << "offset: "; print_simd(offset);
        forward_kernel_simd(source_column_i, offset, current_minimiser, current_offset);
        // std::cout << "current_minimiser: "; print_simd(current_minimiser);
        // std::cout << "current_offset: "; print_simd(current_offset);

        temporary_minimiser_data[i] = current_minimiser;
        temporary_offset_data[i] = current_offset;
    }

    // for (size_t i = 0; i < window_size; ++i)
    // {
    //     print_simd(temporary_minimiser_data[i]);
    // }

    for (size_t i = 0; i < simd_len - 1; ++i)
    {
        // std::cout << "~~~~~~ " << i << std::endl;

        // this might store
        // [x, x, x, x, e, e, e, e]
        //             [x, x, x, x, e, e, e, e]
        //                         [x, x, x, x, e, e, e, e]
        //                                     [x, x, x, x, e, e, e, e]
        int32x8_t target_row_i = simd_gather(simd_len, i, (int const *)&temporary_minimiser_data[0][0]);
        int32x8_t offset_row_i = simd_gather(simd_len, i, (int const *)&temporary_offset_data[0][0]);
        // std::cout << "target_row_i: "; print_simd(target_row_i, window_size);
        // std::cout << "offset_row_i: "; print_simd(offset_row_i, window_size);

        simd_store(target_ptr + (window_size - 1) * i, target_row_i);
        simd_store(offset_ptr + (window_size - 1) * i, offset_row_i);

        assert(((window_size - 1) * (i + 1)) <= target_size);
    }

    for (size_t i = simd_len - 1; i < simd_len; ++i)
    {
        int32x8_t iota{0, 1, 2, 3, 4, 5, 6, 7};
        int32x8_t mask_store{iota < simd_set<int32x8_t>(window_size)};
        int32x8_t target_row_i = simd_gather(simd_len, i, (int const *)&temporary_minimiser_data[0][0]);
        int32x8_t offset_row_i = simd_gather(simd_len, i, (int const *)&temporary_offset_data[0][0]);
        // std::cout << "target_row_i: "; print_simd(target_row_i, window_size);
        // std::cout << "offset_row_i: "; print_simd(offset_row_i, window_size);

        simd_maskstore(target_ptr + (window_size - 1) * i, mask_store, target_row_i);
        simd_maskstore(offset_ptr + (window_size - 1) * i, mask_store, offset_row_i);
    }
}

void compute_backward_simd(size_t const window_size, int const * __restrict source_ptr, int * __restrict target_ptr, size_t const target_size, int * __restrict offset_ptr)
{
    int32x8_t current_minimiser = simd_set<int32x8_t>(std::numeric_limits<int>::max());
    int32x8_t current_offset = simd_set<int32x8_t>(window_size);

    int32x8_t temporary_minimiser_data[8];
    int32x8_t temporary_offset_data[8];

    for (size_t i = 0; i < window_size - 1; ++i)
    {
        // std::cout << "~~~~~~ " << i << std::endl;
        int32x8_t source_column_i = simd_gather(window_size, i, source_ptr);
        // std::cout << "source_column_i: "; print_simd(source_column_i);
        int32x8_t offset = simd_set<int32x8_t>(i);
        // std::cout << "offset: "; print_simd(offset);
        forward_kernel_simd(source_column_i, offset, current_minimiser, current_offset);
        // std::cout << "current_minimiser: "; print_simd(current_minimiser);
        // std::cout << "current_offset: "; print_simd(current_offset);

        temporary_minimiser_data[i] = current_minimiser;
        temporary_offset_data[i] = current_offset;
    }

    // for (size_t i = 0; i < window_size; ++i)
    // {
    //     print_simd(temporary_minimiser_data[i]);
    // }

    for (size_t i = 0; i < simd_len - 1; ++i)
    {
        // std::cout << "~~~~~~ " << i << std::endl;

        // this might store
        // [x, x, x, x, e, e, e, e]
        //             [x, x, x, x, e, e, e, e]
        //                         [x, x, x, x, e, e, e, e]
        //                                     [x, x, x, x, e, e, e, e]
        int32x8_t target_row_i = simd_gather(simd_len, i, (int const *)&temporary_minimiser_data[0][0]);
        int32x8_t offset_row_i = simd_gather(simd_len, i, (int const *)&temporary_offset_data[0][0]);
        // std::cout << "target_row_i: "; print_simd(target_row_i, window_size);
        // std::cout << "offset_row_i: "; print_simd(offset_row_i, window_size);

        simd_store(target_ptr + (window_size - 1) * i, target_row_i);
        simd_store(offset_ptr + (window_size - 1) * i, offset_row_i);

        assert(((window_size - 1) * (i + 1)) <= target_size);
    }

    for (size_t i = simd_len - 1; i < simd_len; ++i)
    {
        int32x8_t iota{0, 1, 2, 3, 4, 5, 6, 7};
        int32x8_t mask_store{iota < simd_set<int32x8_t>(window_size)};
        int32x8_t target_row_i = simd_gather(simd_len, i, (int const *)&temporary_minimiser_data[0][0]);
        int32x8_t offset_row_i = simd_gather(simd_len, i, (int const *)&temporary_offset_data[0][0]);
        // std::cout << "target_row_i: "; print_simd(target_row_i, window_size);
        // std::cout << "offset_row_i: "; print_simd(offset_row_i, window_size);

        simd_maskstore(target_ptr + (window_size - 1) * i, mask_store, target_row_i);
        simd_maskstore(offset_ptr + (window_size - 1) * i, mask_store, offset_row_i);
    }
}

void gather_test()
{
    size_t window_size = 6;
    std::vector<int> values{
        3, 8, 5, 4, 10,17,
        9, 15,20,13,19,13,
        0, 1, 4, 0, 8, 18,
        8, 16,14,9, 14,8,
        19,3, 3, 9, 5, 20,
        0, 1, 2, 3, 4, 5,
        6, 7, 8, 9, 10,11,
        11,10,9, 8, 7, 6};

    {
        int32x8_t column1 = simd_gather(window_size, 0, values.data());
        int column1_expected[]{3, 9, 0, 8, 19, 0, 6, 11};
        assert(std::equal(column1_expected, column1_expected + 8, (int *)(&column1[0])));
    }

    {
        int32x8_t column2 = simd_gather(window_size, 1, values.data());
        int column2_expected[]{8, 15, 1, 16, 3, 1, 7, 10};
        assert(std::equal(column2_expected, column2_expected + 8, (int *)(&column2[0])));
    }

    {
        int32x8_t column3 = simd_gather(window_size, 2, values.data());
        int column3_expected[]{5, 20, 4, 14, 3, 2, 8, 9};
        assert(std::equal(column3_expected, column3_expected + 8, (int *)(&column3[0])));
    }

    {
        int32x8_t column4 = simd_gather(window_size, 3, values.data());
        int column4_expected[]{4, 13, 0, 9, 9, 3, 9, 8};
        assert(std::equal(column4_expected, column4_expected + 8, (int *)(&column4[0])));
    }

    {
        int32x8_t column5 = simd_gather(window_size, 4, values.data());
        int column5_expected[]{10, 19, 8, 14, 5, 4, 10, 7};
        assert(std::equal(column5_expected, column5_expected + 8, (int *)(&column5[0])));
    }

    {
        int32x8_t column6 = simd_gather(window_size, 5, values.data());
        int column6_expected[]{17, 13, 18, 8, 20, 5, 11, 6};
        assert(std::equal(column6_expected, column6_expected + 8, (int *)(&column6[0])));
    }
}

} // namespace stellar::test
