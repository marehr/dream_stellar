
#pragma once

#include "../simd_swizzle_state_rules.hpp"

simd_swizzle_state_rules _mm256_unpacklo_epi32_transitions{
    // ------ _mm256_unpacklo_epi32 ------
    // 0, 1,  2,  3,  4,  5,  6,  7
    // 8, 9, 10, 11, 12, 13, 14, 15
    simd_swizzle_op{{0, 8, 1, 9,/*|*/ 4, 12, 5, 13}, "_mm256_unpacklo_epi32", argument_m256i{}, argument_m256i{}, argument_none{}}//,
};

template <size_t> struct simd_vars;
template <size_t count> struct simd_swizzle_solver;

template <template <size_t> typename simd_swizzle_solver_t = simd_swizzle_solver, template <size_t> typename simd_vars_t = simd_vars>
void simd_swizzle_solver__mm256_unpacklo_epi32_test(auto const & transitions)
{
    // _mm256_unpacklo_epi32
    {
        auto && [memory, final_state] = simd_swizzle_solver_t<2>{}.solve(
            transitions,
            simd_vars_t<2>::iota(), // start state
            simd_vars_t<1>{0, 8, 1, 9,/*|*/ 4, 12, 5, 13}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_m256i, argument_none>>(final_state.ops[0]);
        assert(state0.name == "_mm256_unpacklo_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg1.position == 1);
        assert((memory == simd_vars_t<2>{simd_int32x8{0, 8, 1, 9,/*|*/ 4, 12, 5, 13}, simd_int32x8{8, 9, 10, 11,/*|*/ 12, 13, 14, 15}}));
    }
}
