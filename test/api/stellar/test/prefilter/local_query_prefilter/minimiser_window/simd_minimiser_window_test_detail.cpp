#include <cassert>
#include <iostream>
#include <limits>
#include <span>
#include <vector>

#include <x86intrin.h>

#include <stellar/test/prefilter/local_query_prefilter/minimiser_window/minimiser_state.hpp>

constexpr size_t simd_len = 8;
using int32x8_t [[gnu::vector_size(sizeof(int) * simd_len)]] = int;
using int32x8_mask_t = int32x8_t;

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

inline void simd_store(int * mem_addr, int32x8_t a)
{
    _mm256_storeu_si256((__m256i *)mem_addr, (__m256i)a);
}
inline void simd_maskstore(int * mem_addr, int32x8_t mask,  int32x8_t a)
{
    _mm256_maskstore_epi32((int *)mem_addr, (__m256i)mask, (__m256i)a);
}

inline void forward_kernel_simd(int32x8_t const source, int32x8_t offset, int32x8_t & current_minimiser, int32x8_t & current_offset)
{
    int32x8_mask_t smaller = source < current_minimiser;
    current_minimiser = smaller ? source : current_minimiser;
    current_offset = smaller ? offset : current_offset;
}

void print_simd(int32x8_t simd, size_t limit = simd_len)
{
    std::span<int> values{reinterpret_cast<int(&)[8]>(simd)};
    for (auto & v : values.subspan(0, limit))
        std::cout << v << ", ";
    std::cout << std::endl;
}

void print_chunk(auto * data, size_t const chunk_id, size_t const window_size, size_t const limit)
{
    std::span span{data + chunk_id * window_size, limit};
    for (auto & v : span)
        std::cout << v << ", ";
    std::cout << std::endl;
};

void print_chunked(auto const & data, size_t const window_size)
{
    size_t chunks = data.size() / window_size;
    size_t remaining_elements = data.size() % window_size;
    for (size_t i = 0; i < data.size() / window_size; ++i)
        print_chunk(data.data(), i, window_size, window_size);

    if (remaining_elements > 0)
    {
        print_chunk(data.data(), chunks, window_size, remaining_elements);
    }
};

void compute_forward_simd(size_t const window_size, int const * __restrict source_ptr, int * __restrict target_ptr, int * const __restrict target_end_ptr, int * __restrict offset_ptr)
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

        assert((target_ptr + (window_size - 1) * (i + 1)) <= target_end_ptr);
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

void compute_forward_chunk(size_t const window_size, int const * __restrict source_ptr, int * __restrict target_ptr, int * __restrict offset_ptr)
{
    int current_minimiser = std::numeric_limits<int>::max();
    size_t current_offset = window_size;
    for (size_t offset = window_size; offset > 0;)
    {
        --source_ptr;
        --target_ptr;
        --offset_ptr;
        --offset;

        bool smaller = *source_ptr < current_minimiser;
        current_minimiser = smaller ? *source_ptr : current_minimiser;
        current_offset = smaller ? offset : current_offset;
        *target_ptr = current_minimiser;
        *offset_ptr = current_offset;
    }
}

void compute_backward_chunk(size_t const window_size, int const * __restrict source_ptr, int * __restrict target_ptr, int * __restrict offset_ptr)
{
    int current_minimiser = std::numeric_limits<int>::max();
    size_t current_offset = 0;
    for (size_t offset = 0; offset < window_size; ++offset, ++source_ptr, ++target_ptr, ++offset_ptr)
    {
        bool smaller = *source_ptr < current_minimiser;
        current_minimiser = smaller ? *source_ptr : current_minimiser;
        current_offset = smaller ? offset : current_offset;
        *target_ptr = current_minimiser;
        *offset_ptr = current_offset;
    }
}

void compute_forward_full(size_t const window_size, int const * __restrict source_ptr, int * __restrict target_ptr, size_t const target_size, int * __restrict offset_ptr)
{
    std::ptrdiff_t remaining_elements = target_size % (window_size - 1);
    int * const target_peeled_end_ptr = target_ptr + target_size - remaining_elements;
    for (; true ;)
    {
        source_ptr += window_size - 1;
        target_ptr += window_size - 1;
        offset_ptr += window_size - 1;

        if (target_ptr > target_peeled_end_ptr)
            break;

        compute_forward_chunk(window_size - 1, source_ptr, target_ptr, offset_ptr);
    }

    if (remaining_elements > 0)
    {
        source_ptr -= window_size - 1;
        target_ptr -= window_size - 1;
        offset_ptr -= window_size - 1;
        source_ptr += remaining_elements;
        target_ptr += remaining_elements;
        offset_ptr += remaining_elements;
        compute_forward_chunk(remaining_elements, source_ptr, target_ptr, offset_ptr);
    }
}

void compute_backward_full(size_t const window_size, int const * __restrict source_ptr, int * __restrict target_ptr, size_t const target_size, int * __restrict offset_ptr)
{
    std::ptrdiff_t remaining_elements = target_size % (window_size - 1);
    int * const target_peeled_end_ptr = target_ptr + target_size - remaining_elements;
    for (; target_ptr < target_peeled_end_ptr ;)
    {
        compute_backward_chunk(window_size - 1, source_ptr, target_ptr, offset_ptr);

        source_ptr += window_size - 1;
        target_ptr += window_size - 1;
        offset_ptr += window_size - 1;
    }

    if (remaining_elements > 0)
    {
        // source_ptr -= window_size - 1;
        // target_ptr -= window_size - 1;
        // offset_ptr -= window_size - 1;
        // source_ptr += remaining_elements;
        // target_ptr += remaining_elements;
        // offset_ptr += remaining_elements;
        compute_backward_chunk(remaining_elements, source_ptr, target_ptr, offset_ptr);
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

using stellar::test::minimiser_state;

template <typename value_t>
struct simd_minimiser_window
{
    size_t window_size{};

    using index_t = value_t;

    std::array<value_t, 20> forward_minimizer{};
    std::array<index_t, 20> forward_minimizer_offset{};
    std::array<value_t, 20> backward_minimizer{};
    std::array<index_t, 20> backward_minimizer_offset{};

    struct chunk_state_t
    {
        value_t * forward_it;
        value_t * forward_end;
        value_t * backward_it;
        index_t * forward_offset_it;
        index_t * backward_offset_it;
        value_t minimiser{};
        index_t minimiser_position{};
        bool minimiser_in_forward{false};
    };


    chunk_state_t chunk_state{};

    template <typename iterator_t>
    iterator_t initialize(iterator_t it, iterator_t sentinel)
    {
        std::cout << "values: " << std::endl;

        struct values_t
        {
            iterator_t it;
            iterator_t sentinel;

            size_t size() const
            {
                return sentinel - it;
            }

            value_t const * data() const
            {
                return &*it;
            }
        };
        print_chunked(values_t{it, sentinel}, window_size);

        compute_forward_full(window_size, &*it, forward_minimizer.begin(), forward_minimizer.size(), forward_minimizer_offset.data());
        std::cout << "forward_minimizer: " << std::endl;
        print_chunked(forward_minimizer, window_size - 1);
        std::cout << "forward_minimizer_offset: " << std::endl;
        print_chunked(forward_minimizer_offset, window_size - 1);
        compute_backward_full(window_size, &*it, backward_minimizer.begin(), backward_minimizer.size(), backward_minimizer_offset.data());
        std::cout << "backward_minimizer: " << std::endl;
        print_chunked(backward_minimizer, window_size - 1);
        std::cout << "backward_minimizer_offset: " << std::endl;
        print_chunked(backward_minimizer_offset, window_size - 1);

        chunk_state.forward_it = forward_minimizer.data();
        chunk_state.forward_end = chunk_state.forward_it + window_size - 1;
        chunk_state.backward_it = backward_minimizer.data() + window_size - 1;
        chunk_state.forward_offset_it = forward_minimizer_offset.data();
        chunk_state.backward_offset_it = backward_minimizer_offset.data() + window_size - 1;

        size_t chunk_size = window_size - 1;
        std::cout << "initialize::forward_minimser: ||=" << (chunk_size) << std::endl;
        print_chunked(std::span(chunk_state.forward_it, chunk_size), window_size - 1);
        std::cout << "initialize::backward_minimser: ||=" << (chunk_size) << std::endl;
        print_chunked(std::span(chunk_state.backward_it, chunk_size), window_size - 1);
        chunk_state.minimiser_in_forward = *chunk_state.forward_it < *chunk_state.backward_it;
        chunk_state.minimiser_position = chunk_state.minimiser_in_forward ? *chunk_state.forward_offset_it : window_size - 1 + *chunk_state.backward_offset_it;
        chunk_state.minimiser = std::min(*chunk_state.forward_it, *chunk_state.backward_it);
        return it;
    }

    value_t min() const
    {
        return chunk_state.minimiser;
    }

    minimiser_state cyclic_push()
    {
        ++chunk_state.forward_it;
        ++chunk_state.backward_it;
        ++chunk_state.forward_offset_it;
        ++chunk_state.backward_offset_it;

        std::cout << "forward_minimizer: " << std::endl;
        print_chunked(forward_minimizer, window_size - 1);
        std::cout << "forward_minimizer_offset: " << std::endl;
        print_chunked(forward_minimizer_offset, window_size - 1);
        std::cout << "backward_minimizer: " << std::endl;
        print_chunked(backward_minimizer, window_size - 1);
        std::cout << "backward_minimizer_offset: " << std::endl;
        print_chunked(backward_minimizer_offset, window_size - 1);

        bool new_chunk = chunk_state.forward_it == chunk_state.forward_end;
        if (new_chunk)
        {
            std::cout << "cyclic_push::rechunk: " << std::endl;
            chunk_state.forward_end = chunk_state.forward_it + window_size - 1;
            chunk_state.minimiser_position = chunk_state.minimiser_in_forward ? -1 : chunk_state.minimiser_position - window_size + 1;
        }

        std::ptrdiff_t chunk_size = (chunk_state.forward_end - chunk_state.forward_it);
        assert(chunk_size >= 0);
        std::cout << "cyclic_push::forward_minimser: ||=" << (chunk_size) << std::endl;
        print_chunked(std::span(chunk_state.forward_it, chunk_size), window_size - 1);
        std::cout << "cyclic_push::forward_minimser_offset: ||=" << (chunk_size) << std::endl;
        print_chunked(std::span(chunk_state.forward_offset_it, chunk_size), window_size - 1);
        std::cout << "cyclic_push::backward_minimser: ||=" << (chunk_size) << std::endl;
        print_chunked(std::span(chunk_state.backward_it, chunk_size), window_size - 1);
        std::cout << "cyclic_push::backward_minimser_offset: ||=" << (chunk_size) << std::endl;
        print_chunked(std::span(chunk_state.backward_offset_it, chunk_size), window_size - 1);

        std::ptrdiff_t position = window_size - 1 - chunk_size;
        bool forward_left_window = *(chunk_state.forward_offset_it-1) != *chunk_state.forward_offset_it;
        bool backward_left_window = *(chunk_state.backward_offset_it-1) != *chunk_state.backward_offset_it;
        std::cout << "cyclic_push::position: " << position << std::endl;
        std::cout << "cyclic_push::minimiser_position: " << chunk_state.minimiser_position << std::endl;
        std::cout << "cyclic_push::new_chunk: " << (new_chunk ? "true" : "false") << std::endl;
        std::cout << "cyclic_push::forward_left_window: " << (forward_left_window ? "true" : "false") << std::endl;
        std::cout << "cyclic_push::backward_left_window: " << (backward_left_window ? "true" : "false") << std::endl;

        bool forward_is_smaller_than_backward = *chunk_state.forward_it < *chunk_state.backward_it;
        value_t new_minimizer = forward_is_smaller_than_backward ? *chunk_state.forward_it : *chunk_state.backward_it;
        bool minimiser_changed = new_minimizer != chunk_state.minimiser;
        bool minimiser_position_changed = (position > chunk_state.minimiser_position);
        bool minimiser_left_window = minimiser_position_changed; // && ((chunk_state.minimiser_in_forward && forward_left_window) || (!chunk_state.minimiser_in_forward && backward_left_window));
        chunk_state.minimiser = new_minimizer;
        chunk_state.minimiser_in_forward = true ? forward_is_smaller_than_backward : chunk_state.minimiser_in_forward;
        chunk_state.minimiser_position = forward_is_smaller_than_backward ? *chunk_state.forward_offset_it : window_size - 1 + *chunk_state.backward_offset_it;
        std::cout << "cyclic_push::minimiser_position_changed: " << (minimiser_position_changed ? "true" : "false") << std::endl;
        std::cout << "cyclic_push::minimiser is in : " << (forward_is_smaller_than_backward ? "forward" : "backward") << std::endl;
        std::cout << "cyclic_push::minimiser_changed: " << (minimiser_changed ? "true" : "false") << std::endl;
        std::cout << "cyclic_push::minimiser_left_window: " << (minimiser_left_window ? "true" : "false") << std::endl;
        std::cout << "cyclic_push::minimiser_in_forward: " << (chunk_state.minimiser_in_forward ? "true" : "false") << std::endl;
        std::cout << "cyclic_push::minimiser_position: " << chunk_state.minimiser_position << std::endl;
        std::cout << "cyclic_push::minimiser: " << chunk_state.minimiser << std::endl;

        if (minimiser_left_window)
            return minimiser_state::left_window;
        else if (minimiser_changed)
            return minimiser_state::new_minimizer;
        else
            return minimiser_state::unchanged;
    }
};

void minimiser_test()
{
    std::vector<int> values{
        2, 4, 12, 0,
        15, 13, 10, 9,
        14, 8, 3, 1,
        5, 11, 7, 6};

    simd_minimiser_window<int> minimiser_window{4};

    auto it = minimiser_window.initialize(values.begin(), values.end());


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "[2, 4, 12, *0], 15, 13, 10, 9, 14, 8, 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.min() == 0);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, [4, 12, *0, 15], 13, 10, 9, 14, 8, 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push() == minimiser_state::unchanged);
    assert(minimiser_window.min() == 0);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, [12, *0, 15, 13], 10, 9, 14, 8, 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push() == minimiser_state::unchanged);
    assert(minimiser_window.min() == 0);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, [*0, 15, 13, 10], 9, 14, 8, 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push() == minimiser_state::unchanged);
    assert(minimiser_window.min() == 0);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, [15, 13, 10, *9], 14, 8, 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push() == minimiser_state::left_window);
    assert(minimiser_window.min() == 9);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, [13, 10, *9, 14], 8, 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push() == minimiser_state::unchanged);
    assert(minimiser_window.min() == 9);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, [10, 9, 14, *8], 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push() == minimiser_state::new_minimizer);
    assert(minimiser_window.min() == 8);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, 10, [9, 14, 8, *3], 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push() == minimiser_state::new_minimizer);
    assert(minimiser_window.min() == 3);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, 10, 9, [14, 8, 3, *1], 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push() == minimiser_state::new_minimizer);
    assert(minimiser_window.min() == 1);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, 10, 9, 14, [8, 3, *1, 5], 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push() == minimiser_state::unchanged);
    assert(minimiser_window.min() == 1);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, 10, 9, 14, 8, [3, *1, 5, 11], 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push() == minimiser_state::unchanged);
    assert(minimiser_window.min() == 1);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, 10, 9, 14, 8, 3, [*1, 5, 11, 7], 6" << std::endl;
    assert(minimiser_window.cyclic_push() == minimiser_state::unchanged);
    assert(minimiser_window.min() == 1);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, 10, 9, 14, 8, 3, 1, [*5, 11, 7, 6]" << std::endl;
    assert(minimiser_window.cyclic_push() == minimiser_state::left_window);
    assert(minimiser_window.min() == 5);
}

int main()
{
    minimiser_test();

    std::vector<int> values{
        3 ,8 ,5 ,4 ,10,17,
        9 ,15,20,13,19,13,
        0 ,1 ,4 ,0 ,8 ,18,
        8 ,16,14,9 ,14,8 ,
        19,3 ,3 ,9 ,5 ,20,
        0 ,1 ,2 ,3 ,4 ,5 ,
        6 ,7 ,8 ,9 ,10,11,
        11,10,9 ,8 ,7 ,6 ,
        // 5, 4, 3, 2, 1, 0
    };

    std::vector<int> window_value_forward_expected{
        /*1:*/ 3 /*3 */,4 /*8 */,4 /*5 */,4 /*4 */,10/*10*/,
        /*2:*/ 9 /*17*/,9 /*9 */,13/*15*/,13/*20*/,13/*13*/,
        /*3:*/ 0 /*19*/,0 /*13*/,0 /*0 */,1 /*1 */,4 /*4 */,
        /*4:*/ 0 /*0 */,8 /*8 */,8 /*18*/,8 /*8 */,16/*16*/,
        /*5:*/ 8 /*14*/,8 /*9 */,8 /*14*/,8 /*8 */,19/*19*/,
        /*6:*/ 3 /*3 */,3 /*3 */,5 /*9 */,5 /*5 */,20/*20*/,
        /*7:*/ 0 /*0 */,1 /*1 */,2 /*2 */,3 /*3 */,4 /*4 */,
        /*8:*/ 5 /*5 */,6 /*6 */,7 /*7 */,8 /*8 */,9 /*9 */,
        /*9:*/ 9 /*10*/,9 /*11*/,9 /*11*/,9 /*10*/,9 /*9 */,
        /*0:*/ 6 /*8 */,6 /*7 */,6 /*6 */,
    };
    std::vector<int> window_offset_forward_expected{
        /*1:*/ 0, 3, 3, 3, 4,
        /*2:*/ 1, 1, 4, 4, 4,
        /*3:*/ 2, 2, 2, 3, 4,
        /*4:*/ 0, 3, 3, 3, 4,
        /*5:*/ 3, 3, 3, 3, 4,
        /*6:*/ 1, 1, 3, 3, 4,
        /*7:*/ 0, 1, 2, 3, 4,
        /*8:*/ 0, 1, 2, 3, 4,
        /*9:*/ 4, 4, 4, 4, 4,
        /*0:*/ 2, 2, 2,
    };

    std::vector<int> window_value_backward_expected{
        /*1:*/ 3 /*3 */,3 /*8 */,3 /*5 */,3 /*4 */,3 /*10*/,
        /*2:*/ 17/*17*/,9 /*9 */,9 /*15*/,9 /*20*/,9 /*13*/,
        /*3:*/ 19/*19*/,13/*13*/,0 /*0 */,0 /*1 */,0 /*4 */,
        /*4:*/ 0 /*0 */,0 /*8 */,0 /*18*/,0 /*8 */,0 /*16*/,
        /*5:*/ 14/*14*/,9 /*9 */,9 /*14*/,8 /*8 */,8 /*19*/,
        /*6:*/ 3 /*3 */,3 /*3 */,3 /*9 */,3 /*5 */,3 /*20*/,
        /*7:*/ 0 /*0 */,0 /*1 */,0 /*2 */,0 /*3 */,0 /*4 */,
        /*8:*/ 5 /*5 */,5 /*6 */,5 /*7 */,5 /*8 */,5 /*9 */,
        /*9:*/ 10/*10*/,10/*11*/,10/*11*/,10/*10*/,9 /*9 */,
        /*0:*/ 8 /*8 */,7 /*7 */,6 /*6 */,
    };
    std::vector<int> window_offset_backward_expected{
        /*1:*/ 0, 0, 0, 0, 0,
        /*2:*/ 0, 1, 1, 1, 1,
        /*3:*/ 0, 1, 2, 2, 2,
        /*4:*/ 0, 0, 0, 0, 0,
        /*5:*/ 0, 1, 1, 3, 3,
        /*6:*/ 0, 0, 0, 0, 0,
        /*7:*/ 0, 0, 0, 0, 0,
        /*8:*/ 0, 0, 0, 0, 0,
        /*9:*/ 0, 0, 0, 0, 4,
        /*0:*/ 0, 1, 2,
        // 5, 5, 5, 5, 5, 5,
    };

    std::size_t window_size = 6;

    gather_test();

    std::vector<int> window_value_forward;
    std::vector<int> window_offset_forward;
    window_value_forward.resize(values.size());
    window_offset_forward.resize(values.size());

    std::vector<int> window_value_backward;
    std::vector<int> window_offset_backward;
    window_value_backward.resize(values.size());
    window_offset_backward.resize(values.size());

    compute_forward_full(window_size, values.data(), window_value_forward.data(), window_value_forward.size(), window_offset_forward.data());
    compute_backward_full(window_size, values.data(), window_value_backward.data(), window_value_backward.size(), window_offset_backward.data());

    std::cout << "window_value_forward: " << std::endl;
    print_chunked(window_value_forward, window_size - 1);
    assert(std::equal(window_value_forward_expected.begin(), window_value_forward_expected.end(), window_value_forward.begin()));

    std::cout << "window_offset_forward: " << std::endl;
    print_chunked(window_offset_forward, window_size - 1);
    assert(std::equal(window_offset_forward_expected.begin(), window_offset_forward_expected.end(), window_offset_forward.begin()));

    std::cout << "window_value_backward: " << std::endl;
    print_chunked(window_value_backward, window_size - 1);
    assert(std::equal(window_value_backward_expected.begin(), window_value_backward_expected.end(), window_value_backward.begin()));

    std::cout << "window_offset_backward: " << std::endl;
    print_chunked(window_offset_backward, window_size - 1);
    assert(std::equal(window_offset_backward_expected.begin(), window_offset_backward_expected.end(), window_offset_backward.begin()));

}
