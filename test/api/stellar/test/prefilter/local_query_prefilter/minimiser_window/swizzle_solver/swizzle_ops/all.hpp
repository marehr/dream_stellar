
#pragma once

#include "permute2f128_si256_transitions.hpp"

simd_swizzle_state_rules all_transitions{[]()
{
    simd_swizzle_state_rules _transitions{};
    auto append = [&](auto && other)
    {
        _transitions.transitions.insert(
            _transitions.transitions.end(),
            other.transitions.begin(),
            other.transitions.end()
        );
    };

    append(permute2f128_si256_transitions);

    return _transitions;
}()};

template <size_t> struct simd_vars;
template <size_t count> struct simd_swizzle_solver;

template <template <size_t> typename simd_swizzle_solver_t = simd_swizzle_solver, template <size_t> typename simd_vars_t = simd_vars>
void simd_swizzle_solver_identity_test(auto const & transitions)
{
    // identity
    {
        auto && [memory, final_state] = simd_swizzle_solver_t<1>{}.find(
            transitions,
            simd_vars_t<1>{0, 1, 2, 3, 4, 5, 6, 7}, // start state
            simd_vars_t<1>{0, 1, 2, 3, 4, 5, 6, 7}, // end state
            1
        );

        assert(final_state.size() == 0);
        assert((memory == simd_vars_t<1>{0, 1, 2, 3, 4, 5, 6, 7}));
    }
}
