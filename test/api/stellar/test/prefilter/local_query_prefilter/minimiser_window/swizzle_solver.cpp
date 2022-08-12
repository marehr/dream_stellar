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
}
