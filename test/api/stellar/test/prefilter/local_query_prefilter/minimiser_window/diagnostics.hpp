#pragma once

#include <cassert>
#include <iostream>
#include <span>

namespace stellar::test
{

template <typename simd_t>
void print_simd(simd_t simd, size_t limit)
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

} // namespace stellar::test
