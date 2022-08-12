
#pragma once

#include <array>
#include <iostream>

void print_simd(std::basic_ostream<char> & cout, std::array<int, 8> const & simd)
{
    cout << "[";
    for (int i = 0; i < 8; ++i)
    {
        if (i != 0) cout << ",";
        cout << simd[i];
    }
    cout << "]";
}

struct simd_int32x8 : public std::array<int, 8>
{
    template <typename char_t>
    friend std::basic_ostream<char_t> & operator<<(std::basic_ostream<char_t> & cout, simd_int32x8 const & simd)
    {
        print_simd(cout, simd);
        return cout;
    }
};
