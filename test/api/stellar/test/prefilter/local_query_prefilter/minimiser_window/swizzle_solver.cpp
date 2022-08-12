#include "swizzle_solver/simd_int32x.hpp"
#include "swizzle_solver/simd_swizzle_op.hpp"
#include "swizzle_solver/simd_swizzle_state_rules.hpp"
#include "swizzle_solver/simd_swizzle_op_trace.hpp"
#include "swizzle_solver/simd_vars.hpp"
#include "swizzle_solver/simd_swizzle_solver.hpp"
#include "swizzle_solver/swizzle_ops/all.hpp"

int main()
{
    simd_vars_test();
    simd_swizzle_solver_identity_test(all_transitions);
    simd_swizzle_solver_permute2f128_si256_test(all_transitions);
    simd_swizzle_solver_shuffle_epi32_test(all_transitions);
    simd_swizzle_solver__mm256_unpackhi_epi32_test(all_transitions);
    simd_swizzle_solver__mm256_unpacklo_epi32_test(all_transitions);

    //  0,  1,  2,  3
    //  4,  5,  6,  7
    //  8,  9, 10, 11
    // 12, 13, 14, 15
    // ->
    //  0,  4,  8, 12,
    //  1,  5,  9, 13,
    //  2,  6, 10, 14,
    //  3,  7, 11, 15,
    simd_swizzle_solver<3> solver{};
    auto start_memory = simd_vars<3>::iota();
    simd_vars<2> end_data{simd_int32x8{0, 4, 8, 12, 1, 5, 9, 13}, simd_int32x8{2, 6, 10, 14, 3, 7, 11, 15}};
    auto && [end_memory, traces] = solver.solve(all_transitions, start_memory, end_data, 5);

    std::cout << end_memory << "; " << traces << std::endl;
}
