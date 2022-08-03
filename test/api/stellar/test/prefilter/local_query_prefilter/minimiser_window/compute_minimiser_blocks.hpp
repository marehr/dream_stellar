
#pragma once

#include <cmath>
#include <limits>

namespace stellar::test
{

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
        bool smaller = *source_ptr <= current_minimiser;
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

} // namespace stellar::test
