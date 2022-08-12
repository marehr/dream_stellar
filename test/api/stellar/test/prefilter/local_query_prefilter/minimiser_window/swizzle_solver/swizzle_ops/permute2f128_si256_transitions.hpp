
#pragma once

#include "../simd_swizzle_state_rules.hpp"

simd_swizzle_state_rules permute2f128_si256_transitions{
    // ------ permute2f128_si256 ------
    /*permute2f128_si256_01: A[0], A[0]*/ simd_swizzle_op{{0, 1, 2, 3,/*|*/ 0, 1, 2, 3}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_any_m256i{}, argument_imm8{0x00}},
    // /*permute2f128_si256_02: A[0], A[1]*/ simd_swizzle_op{{0, 1, 2, 3,/*|*/ 4, 5, 6, 7}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_any_m256i{}, argument_imm8{0x10}},
    /*permute2f128_si256_03: A[0], B[0]*/ simd_swizzle_op{{0, 1, 2, 3,/*|*/ 8, 9, 10, 11}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x20}},
    /*permute2f128_si256_04: A[0], B[1]*/ simd_swizzle_op{{0, 1, 2, 3,/*|*/ 12, 13, 14, 15}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x30}},
    /*permute2f128_si256_05: A[1], A[0]*/ simd_swizzle_op{{4, 5, 6, 7,/*|*/ 0, 1, 2, 3}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_any_m256i{}, argument_imm8{0x01}},
    /*permute2f128_si256_06: A[1], A[1]*/ simd_swizzle_op{{4, 5, 6, 7,/*|*/ 4, 5, 6, 7}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_any_m256i{}, argument_imm8{0x11}},
    /*permute2f128_si256_07: A[1], B[0]*/ simd_swizzle_op{{4, 5, 6, 7,/*|*/ 8, 9, 10, 11}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x21}},
    /*permute2f128_si256_08: A[1], B[1]*/ simd_swizzle_op{{4, 5, 6, 7,/*|*/ 12, 13, 14, 15}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x31}},

    // /*permute2f128_si256_09: B[0], A[0]*/ simd_swizzle_op{{8, 9, 10, 11,/*|*/ 0, 1, 2, 3}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x02}},
    // /*permute2f128_si256_10: B[0], A[1]*/ simd_swizzle_op{{8, 9, 10, 11,/*|*/ 4, 5, 6, 7}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x12}},
    // /*permute2f128_si256_11: B[0], B[0]*/ simd_swizzle_op{{8, 9, 10, 11,/*|*/ 8, 9, 10, 11}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x22}},
    // /*permute2f128_si256_12: B[0], B[1]*/ simd_swizzle_op{{8, 9, 10, 11,/*|*/ 12, 13, 14, 15}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x32}},
    // /*permute2f128_si256_13: B[1], A[0]*/ simd_swizzle_op{{12, 13, 14, 15,/*|*/ 0, 1, 2, 3}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x03}},
    // /*permute2f128_si256_14: B[1], A[1]*/ simd_swizzle_op{{12, 13, 14, 15,/*|*/ 4, 5, 6, 7}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x13}},
    // /*permute2f128_si256_15: B[1], B[0]*/ simd_swizzle_op{{12, 13, 14, 15,/*|*/ 8, 9, 10, 11}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x23}},
    // /*permute2f128_si256_16: B[1], B[1]*/ simd_swizzle_op{{12, 13, 14, 15,/*|*/ 12, 13, 14, 15}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x33}},
};

template <size_t> struct simd_vars;
template <size_t count> struct simd_swizzle_solver;

template <template <size_t> typename simd_swizzle_solver_t = simd_swizzle_solver, template <size_t> typename simd_vars_t = simd_vars>
void simd_swizzle_solver_permute2f128_si256_test(auto const & transitions)
{
    // permute2f128_si256_0x00
    {
        auto && [memory, final_state] = simd_swizzle_solver_t<1>{}.solve(
            transitions,
            simd_vars_t<1>::iota(), // start state
            simd_vars_t<1>{0, 1, 2, 3, 0, 1, 2, 3}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_any_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg2.constant == 0x00);
        assert((memory == simd_vars_t<1>{0, 1, 2, 3, 0, 1, 2, 3}));
    }

    // permute2f128_si256_0x01
    {
        auto && [memory, final_state] = simd_swizzle_solver_t<1>{}.solve(
            transitions,
            simd_vars_t<1>::iota(), // start state
            simd_vars_t<1>{4, 5, 6, 7, 0, 1, 2, 3}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_any_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg2.constant == 0x01);
        assert((memory == simd_vars_t<1>{4, 5, 6, 7, 0, 1, 2, 3}));
    }

    // permute2f128_si256_0x11
    {
        auto && [memory, final_state] = simd_swizzle_solver_t<1>{}.solve(
            transitions,
            simd_vars_t<1>::iota(), // start state
            simd_vars_t<1>{4, 5, 6, 7, 4, 5, 6, 7}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_any_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg2.constant == 0x11);
        assert((memory == simd_vars_t<1>{4, 5, 6, 7, 4, 5, 6, 7}));
    }

    // permute2f128_si256_0x20
    {
        std::cout << "~~~~~~~~~~" << std::endl;
        auto && [memory, final_state] = simd_swizzle_solver_t<2>{}.solve(
            transitions,
            simd_vars_t<2>::iota(), // start state
            simd_vars_t<1>{0, 1, 2, 3, 8, 9, 10, 11}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg1.position == 1);
        assert(state0.arg2.constant == 0x20);
        assert((memory == simd_vars_t<2>{simd_int32x8{0,1,2,3,8,9,10,11}, simd_int32x8{8,9,10,11,12,13,14,15}}));
    }

    // permute2f128_si256_0x30
    {
        std::cout << "~~~~~~~~~~" << std::endl;
        auto && [memory, final_state] = simd_swizzle_solver_t<2>{}.solve(
            transitions,
            simd_vars_t<2>::iota(), // start state
            simd_vars_t<1>{0, 1, 2, 3, 12, 13, 14, 15}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg1.position == 1);
        assert(state0.arg2.constant == 0x30);
        assert((memory == simd_vars_t<2>{simd_int32x8{0, 1, 2, 3, 12, 13, 14, 15}, simd_int32x8{8,9,10,11,12,13,14,15}}));
    }

    // permute2f128_si256_0x31
    {
        std::cout << "~~~~~~~~~~" << std::endl;
        auto && [memory, final_state] = simd_swizzle_solver_t<2>{}.solve(
            transitions,
            simd_vars_t<2>::iota(), // start state
            simd_vars_t<1>{4, 5, 6, 7, 12, 13, 14, 15}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg1.position == 1);
        assert(state0.arg2.constant == 0x31);
        assert((memory == simd_vars_t<2>{simd_int32x8{4, 5, 6, 7, 12, 13, 14, 15}, simd_int32x8{8,9,10,11,12,13,14,15}}));
    }

    // permute2f128_si256_0x20 reverse
    {
        std::cout << "~~~~~~~~~~" << std::endl;
        auto && [memory, final_state] = simd_swizzle_solver_t<2>{}.solve(
            transitions,
            simd_vars_t<2>::iota(), // start state
            simd_vars_t<1>{8, 9, 10, 11, 0, 1, 2, 3}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 1);
        assert(state0.arg1.position == 0);
        assert(state0.arg2.constant == 0x20);
        assert((memory == simd_vars_t<2>{simd_int32x8{8, 9, 10, 11, 0, 1, 2, 3}, simd_int32x8{8,9,10,11,12,13,14,15}}));
    }

    // permute2f128_si256_0x21 reverse
    {
        std::cout << "~~~~~~~~~~" << std::endl;
        auto && [memory, final_state] = simd_swizzle_solver_t<2>{}.solve(
            transitions,
            simd_vars_t<2>::iota(), // start state
            simd_vars_t<1>{12, 13, 14, 15, 0, 1, 2, 3}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 1);
        assert(state0.arg1.position == 0);
        assert(state0.arg2.constant == 0x21);
        assert((memory == simd_vars_t<2>{simd_int32x8{12, 13, 14, 15, 0, 1, 2, 3}, simd_int32x8{8,9,10,11,12,13,14,15}}));
    }
}
