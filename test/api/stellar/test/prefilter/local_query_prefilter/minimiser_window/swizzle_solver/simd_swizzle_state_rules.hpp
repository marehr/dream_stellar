
#pragma once

#include <vector>

#include "simd_swizzle_op.hpp"

struct simd_swizzle_state_rules
{
    simd_swizzle_state_rules() = default;

    simd_swizzle_state_rules(simd_swizzle_state_rules const &) = default;

    template <typename ...args_t>
        requires (sizeof...(args_t) > 1)
    simd_swizzle_state_rules(args_t && ...args)
        : transitions{simd_swizzle_opv{std::forward<args_t>(args)}...}
    {}

    std::vector<simd_swizzle_opv> transitions;
};
