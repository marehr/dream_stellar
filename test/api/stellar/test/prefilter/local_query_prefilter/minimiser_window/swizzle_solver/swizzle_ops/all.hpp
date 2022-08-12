
#pragma once

#include "permute2f128_si256_transitions.hpp"
#include "shuffle_epi32_transitions.hpp"
#include "_mm256_unpackhi_epi32.hpp"
#include "_mm256_unpacklo_epi32.hpp"

// YES:
// __m256i _mm256_blend_epi32 (__m256i a, __m256i b, const int imm8) vpblendd
// __m256i _mm256_permute2f128_si256 (__m256i a, __m256i b, int imm8) vperm2f128
// __m256i _mm256_permutevar8x32_epi32 (__m256i a, __m256i idx) vpermd
// __m256i _mm256_shuffle_epi32 (__m256i a, const int imm8) vpshufd
// __m256i _mm256_unpackhi_epi32 (__m256i a, __m256i b) vpunpckhdq
// __m256i _mm256_unpacklo_epi32 (__m256i a, __m256i b) vpunpckldq
//
// MAYBE:
// __m256i _mm256_blendv_epi8 (__m256i a, __m256i b, __m256i mask) vpblendvb
// __m256i _mm256_permute2x128_si256 (__m256i a, __m256i b, const int imm8) vperm2i128
//
// NO (32BIT, but float)
// __m256 _mm256_blend_ps (__m256 a, __m256 b, const int imm8) vblendps
// __m256 _mm256_blendv_ps (__m256 a, __m256 b, __m256 mask) vblendvps
// __m256 _mm256_permute_ps (__m256 a, int imm8) vpermilps
// __m256 _mm256_permute2f128_ps (__m256 a, __m256 b, int imm8) vperm2f128
// __m256 _mm256_permutevar_ps (__m256 a, __m256i b) vpermilps
// __m256 _mm256_permutevar8x32_ps (__m256 a, __m256i idx) vpermps
// __m256 _mm256_shuffle_ps (__m256 a, __m256 b, const int imm8) vshufps
// __m256 _mm256_unpackhi_ps (__m256 a, __m256 b) vunpckhps
// __m256 _mm256_unpacklo_ps (__m256 a, __m256 b) vunpcklps
//
// NO (NOT32BIT):
// __m256i _mm256_blend_epi16 (__m256i a, __m256i b, const int imm8)
// __m256d _mm256_blend_pd (__m256d a, __m256d b, const int imm8) vblendpd
// __m256d _mm256_blendv_pd (__m256d a, __m256d b, __m256d mask) vblendvpd
//
// __m256d _mm256_broadcast_pd (__m128d const * mem_addr) vbroadcastf128
// __m256 _mm256_broadcast_ps (__m128 const * mem_addr) vbroadcastf128
// __m256d _mm256_broadcast_sd (double const * mem_addr) vbroadcastsd
//
// __m256i _mm256_broadcastb_epi8 (__m128i a) vpbroadcastb
// __m256i _mm256_broadcastd_epi32 (__m128i a) vpbroadcastd
// __m256i _mm256_broadcastq_epi64 (__m128i a) vpbroadcastq
// __m256d _mm256_broadcastsd_pd (__m128d a) vbroadcastsd
//
// __m256i _mm_broadcastsi128_si256 (__m128i a) vbroadcasti128
// __m256i _mm256_broadcastsi128_si256 (__m128i a) vbroadcasti128
// __m256 _mm256_broadcastss_ps (__m128 a) vbroadcastss
// __m256i _mm256_broadcastw_epi16 (__m128i a) vpbroadcastw
//
// int _mm256_extract_epi16 (__m256i a, const int index) ...
// __int32 _mm256_extract_epi32 (__m256i a, const int index) ...
// __int64 _mm256_extract_epi64 (__m256i a, const int index) ...
// int _mm256_extract_epi8 (__m256i a, const int index) ...
// __m128d _mm256_extractf128_pd (__m256d a, const int imm8) vextractf128
// __m128 _mm256_extractf128_ps (__m256 a, const int imm8) vextractf128
// __m128i _mm256_extractf128_si256 (__m256i a, const int imm8) vextractf128
// __m128i _mm256_extracti128_si256 (__m256i a, const int imm8) vextracti128
//
// __m256i _mm256_insert_epi16 (__m256i a, __int16 i, const int index) ...
// __m256i _mm256_insert_epi32 (__m256i a, __int32 i, const int index) ...
// __m256i _mm256_insert_epi64 (__m256i a, __int64 i, const int index) ...
// __m256i _mm256_insert_epi8 (__m256i a, __int8 i, const int index) ...
// __m256d _mm256_insertf128_pd (__m256d a, __m128d b, int imm8) vinsertf128
// __m256 _mm256_insertf128_ps (__m256 a, __m128 b, int imm8) vinsertf128
// __m256i _mm256_insertf128_si256 (__m256i a, __m128i b, int imm8) vinsertf128
// __m256i _mm256_inserti128_si256 (__m256i a, __m128i b, const int imm8) vinserti128
//
// __m256d _mm256_permute_pd (__m256d a, int imm8) vpermilpd
// __m256d _mm256_permute2f128_pd (__m256d a, __m256d b, int imm8) vperm2f128
//
// __m256i _mm256_permute4x64_epi64 (__m256i a, const int imm8) vpermq
// __m256d _mm256_permute4x64_pd (__m256d a, const int imm8) vpermpd
// __m256d _mm256_permutevar_pd (__m256d a, __m256i b) vpermilpd
//
// __m256i _mm256_shuffle_epi8 (__m256i a, __m256i b) vpshufb
// __m256d _mm256_shuffle_pd (__m256d a, __m256d b, const int imm8) vshufpd
//
// __m256i _mm256_shufflehi_epi16 (__m256i a, const int imm8) vpshufhw
// __m256i _mm256_shufflelo_epi16 (__m256i a, const int imm8) vpshuflw
// __m256i _mm256_unpackhi_epi16 (__m256i a, __m256i b) vpunpckhwd
//
// __m256i _mm256_unpackhi_epi64 (__m256i a, __m256i b) vpunpckhqdq
// __m256i _mm256_unpackhi_epi8 (__m256i a, __m256i b) vpunpckhbw
// __m256d _mm256_unpackhi_pd (__m256d a, __m256d b) vunpckhpd
// __m256i _mm256_unpacklo_epi16 (__m256i a, __m256i b) vpunpcklwd
// __m256i _mm256_unpacklo_epi64 (__m256i a, __m256i b) vpunpcklqdq
// __m256i _mm256_unpacklo_epi8 (__m256i a, __m256i b) vpunpcklbw
// __m256d _mm256_unpacklo_pd (__m256d a, __m256d b) vunpcklpd

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

    append(_mm256_unpackhi_epi32_transitions);
    append(_mm256_unpacklo_epi32_transitions);
    append(permute2f128_si256_transitions);
    append(shuffle_epi32_transitions);

    return _transitions;
}()};

template <size_t> struct simd_vars;
template <size_t count> struct simd_swizzle_solver;

template <template <size_t> typename simd_swizzle_solver_t = simd_swizzle_solver, template <size_t> typename simd_vars_t = simd_vars>
void simd_swizzle_solver_identity_test(auto const & transitions)
{
    // identity
    {
        auto && [memory, final_state] = simd_swizzle_solver_t<1>{}.solve(
            transitions,
            simd_vars_t<1>{0, 1, 2, 3, 4, 5, 6, 7}, // start state
            simd_vars_t<1>{0, 1, 2, 3, 4, 5, 6, 7}, // end state
            1
        );

        assert(final_state.size() == 0);
        assert((memory == simd_vars_t<1>{0, 1, 2, 3, 4, 5, 6, 7}));
    }
}
