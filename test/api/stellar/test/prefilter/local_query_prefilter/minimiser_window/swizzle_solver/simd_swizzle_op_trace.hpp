
#pragma once

#include <vector>

#include "simd_swizzle_op.hpp"

struct simd_swizzle_op_trace
{
    size_t size() const
    {
        return ops.size();
    }

    std::vector<simd_swizzle_opv> ops;

    template <typename char_t>
    friend std::basic_ostream<char_t> & operator<<(std::basic_ostream<char_t> & cout, simd_swizzle_op_trace const & trace)
    {
        cout << "[";
        size_t i = 0;
        for (simd_swizzle_opv const & opv : trace.ops)
        {
            std::visit([&cout](auto const & op)
            {
                cout << op;
            }, opv);

            if (++i < trace.ops.size())
                cout << " -> ";
        }
        cout << "]";
        return cout;
    }
};
