
#pragma once

#include <cmath>
#include <vector>

#include <stellar/test/prefilter/local_query_prefilter/minimiser_window/minimiser_state.hpp>

#include "compute_minimiser_blocks.hpp"
#include "diagnostics.hpp"

namespace stellar::test
{

template <typename value_t>
struct simd_minimiser_window
{
    size_t window_size{};

    using index_t = value_t;

    std::vector<value_t> _memory{};

    std::span<value_t> forward_minimizer{};
    std::span<index_t> forward_minimizer_offset{};
    std::span<value_t> backward_minimizer{};
    std::span<index_t> backward_minimizer_offset{};

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
        iterator_t begin = it;
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

        size_t size = (sentinel - it) + 1;
        _memory.reserve(4 * size);
        forward_minimizer = {_memory.data(), size};
        forward_minimizer_offset = {_memory.data() + size, size};
        backward_minimizer = {_memory.data() + 2 * size, size};
        backward_minimizer_offset = {_memory.data() + 3 * size, size};

        assert(sentinel - it < forward_minimizer.size());

        compute_forward_full(window_size, &*it, forward_minimizer.data(), forward_minimizer.size(), forward_minimizer_offset.data());
        std::cout << "forward_minimizer: " << std::endl;
        print_chunked(forward_minimizer, window_size - 1);
        std::cout << "forward_minimizer_offset: " << std::endl;
        print_chunked(forward_minimizer_offset, window_size - 1);
        compute_backward_full(window_size, &*it, backward_minimizer.data(), backward_minimizer.size(), backward_minimizer_offset.data());
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
        return begin + window_size - 1;
    }

    value_t min() const
    {
        return chunk_state.minimiser;
    }

    minimiser_state cyclic_push(value_t const new_value)
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
        bool minimiser_position_changed = (position > chunk_state.minimiser_position);
        bool minimiser_left_window = minimiser_position_changed; // && ((chunk_state.minimiser_in_forward && forward_left_window) || (!chunk_state.minimiser_in_forward && backward_left_window));
        bool minimiser_changed = minimiser_left_window || new_minimizer != chunk_state.minimiser;
        chunk_state.minimiser = new_minimizer;
        chunk_state.minimiser_in_forward = minimiser_changed ? forward_is_smaller_than_backward : chunk_state.minimiser_in_forward;
        index_t new_minimizer_position = forward_is_smaller_than_backward ? *chunk_state.forward_offset_it : window_size - 1 + *chunk_state.backward_offset_it;
        chunk_state.minimiser_position = minimiser_changed ? new_minimizer_position : chunk_state.minimiser_position;
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

using stellar::test::minimiser_state;

void minimiser_test()
{
    int values[] = {
        2, 4, 12, 0,
        15, 13, 10, 9,
        14, 8, 3, 1,
        5, 11, 7, 6};

    simd_minimiser_window<int> minimiser_window{4};

    auto it = minimiser_window.initialize(std::begin(values), std::end(values));

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "[2, 4, 12, *0], 15, 13, 10, 9, 14, 8, 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.min() == 0);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, [4, 12, *0, 15], 13, 10, 9, 14, 8, 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push(0) == minimiser_state::unchanged);
    assert(minimiser_window.min() == 0);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, [12, *0, 15, 13], 10, 9, 14, 8, 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push(0) == minimiser_state::unchanged);
    assert(minimiser_window.min() == 0);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, [*0, 15, 13, 10], 9, 14, 8, 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push(0) == minimiser_state::unchanged);
    assert(minimiser_window.min() == 0);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, [15, 13, 10, *9], 14, 8, 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push(0) == minimiser_state::left_window);
    assert(minimiser_window.min() == 9);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, [13, 10, *9, 14], 8, 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push(0) == minimiser_state::unchanged);
    assert(minimiser_window.min() == 9);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, [10, 9, 14, *8], 3, 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push(0) == minimiser_state::new_minimizer);
    assert(minimiser_window.min() == 8);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, 10, [9, 14, 8, *3], 1, 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push(0) == minimiser_state::new_minimizer);
    assert(minimiser_window.min() == 3);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, 10, 9, [14, 8, 3, *1], 5, 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push(0) == minimiser_state::new_minimizer);
    assert(minimiser_window.min() == 1);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, 10, 9, 14, [8, 3, *1, 5], 11, 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push(0) == minimiser_state::unchanged);
    assert(minimiser_window.min() == 1);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, 10, 9, 14, 8, [3, *1, 5, 11], 7, 6" << std::endl;
    assert(minimiser_window.cyclic_push(0) == minimiser_state::unchanged);
    assert(minimiser_window.min() == 1);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, 10, 9, 14, 8, 3, [*1, 5, 11, 7], 6" << std::endl;
    assert(minimiser_window.cyclic_push(0) == minimiser_state::unchanged);
    assert(minimiser_window.min() == 1);


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "2, 4, 12, 0, 15, 13, 10, 9, 14, 8, 3, 1, [*5, 11, 7, 6]" << std::endl;
    assert(minimiser_window.cyclic_push(0) == minimiser_state::left_window);
    assert(minimiser_window.min() == 5);
}

} // namespace stellar::test
