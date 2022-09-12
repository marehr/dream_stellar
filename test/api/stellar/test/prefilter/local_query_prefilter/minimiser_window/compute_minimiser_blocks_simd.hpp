
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

} // namespace stellar::test
