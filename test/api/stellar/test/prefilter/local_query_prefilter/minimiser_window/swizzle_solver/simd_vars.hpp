
#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <unordered_map>

#include "simd_int32x.hpp"

template <size_t count>
struct simd_vars : public std::array<simd_int32x8, count>
{
    static constexpr simd_vars iota()
    {
        simd_vars _iota{};

        for (size_t offset = 0, value = 0; offset < count; ++offset)
        {
            for (size_t i = 0; i < _iota[offset].size(); ++i, ++value)
            {
                _iota[offset][i] = value;
                // std::cout << "_iota[" << offset << "][" << i << "] = i + offset * count = " << (i + offset * count) << std::endl;
            }
        }

        return _iota;
    }

    struct hasher
    {
        std::size_t operator()(simd_vars const& var) const noexcept
        {
            size_t hash{};
            for (simd_int32x8 memory : var)
            {
                for (int idx : memory)
                {
                    hash = hash ^ (std::hash<int>{}(idx) << 1);
                }
            }
            return hash; // or use boost::hash_combine
        }
    };

    template <size_t count2>
    bool contains(simd_vars<count2> const & memory) const
    {
        simd_vars<count> const & memory1{*this};
        simd_vars<count2> const & memory2{memory};

        for (simd_int32x8 const & simd : memory2)
        {
            auto it = std::find(memory1.begin(), memory1.end(), simd);
            if (it == memory1.end())
                return false;
        }

        return true;
    }

    template <size_t count2>
    bool has_at_least_each_once(simd_vars<count2> const & final_memory) const
    {
        simd_vars<count> const & memory1{*this};
        simd_vars<count2> const & memory2{final_memory};

        std::unordered_map<int, int> map{};

        for (simd_int32x8 const & simd : memory1)
        {
            for (int v : simd)
                ++map[v];
        }

        for (simd_int32x8 const & simd : memory2)
        {
            for (int v : simd)
                if (!map.contains(v))
                    return false;
        }

        return true;
    }

    template <typename char_t>
    friend std::basic_ostream<char_t> & operator<<(std::basic_ostream<char_t> & cout, simd_vars<count> const & vars)
    {
        for (size_t i = 0; i < count; ++i)
        {
            simd_int32x8 const & simd = vars[i];
            cout << simd;
            if (i < count - 1)
                cout << std::endl;
        }
        return cout;
    }
};


void simd_vars_test()
{
    assert((simd_vars<1>::iota() == simd_vars<1>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}}));
    assert((simd_vars<2>::iota() == simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{8+0, 8+1, 8+2, 8+3, 8+4, 8+5, 8+6, 8+7}}));
    assert((simd_vars<3>::iota() == simd_vars<3>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{8+0, 8+1, 8+2, 8+3, 8+4, 8+5, 8+6, 8+7}, simd_int32x8{16+0, 16+1, 16+2, 16+3, 16+4, 16+5, 16+6, 16+7}}));

    assert(!(simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}}.contains(simd_vars<1>{{0, 1, 2, 3, 4, 5, 7, 6}})));
    assert((simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}}.contains(simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}})));
    assert((simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{-1, -1, -1, -1, -1, -1, -1, -1}}.contains(simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}})));
    assert((simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{8, 9, 10, 11, 12, 13, 14, 15}}.contains(simd_vars<1>{{8, 9, 10, 11, 12, 13, 14, 15}})));
    assert(!(simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{8, 9, 10, 11, 12, 13, 14, 15}}.contains(simd_vars<1>{{4, 5, 6, 7, 8, 9, 10, 11}})));

    assert((simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}}.has_at_least_each_once(simd_vars<1>{{0, 1, 2, 3, 4, 5, 7, 6}})));
    assert((simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}}.has_at_least_each_once(simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}})));
    assert((simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{-1, -1, -1, -1, -1, -1, -1, -1}}.has_at_least_each_once(simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}})));
    assert((simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{8, 9, 10, 11, 12, 13, 14, 15}}.has_at_least_each_once(simd_vars<1>{{8, 9, 10, 11, 12, 13, 14, 15}})));
    assert((simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{8, 9, 10, 11, 12, 13, 14, 15}}.has_at_least_each_once(simd_vars<1>{{4, 5, 6, 7, 8, 9, 10, 11}})));
    assert(!(simd_vars<2>{simd_int32x8{0, 1, 2, 3, 5, 5, 6, 7}, simd_int32x8{8, 9, 10, 11, 12, 13, 14, 15}}.has_at_least_each_once(simd_vars<1>{{4, 5, 6, 7, 8, 9, 10, 11}})));
}
