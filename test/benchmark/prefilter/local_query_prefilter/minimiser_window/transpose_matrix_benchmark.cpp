#include <benchmark/benchmark.h>

#include <random>

#include <seqan3/utility/container/aligned_allocator.hpp>

#include "../test/api/stellar/test/prefilter/local_query_prefilter/minimiser_window/transpose_matrix.hpp"

using value_type = int;
static constexpr size_t values_size = 100'000;

std::vector<value_type> _g_values{[]()
{
    std::random_device random_device{};
    std::size_t random_seed = random_device();
    std::mt19937_64 random_engine(random_seed);
    std::uniform_int_distribution<value_type> value_distribution(0, 20);

    std::vector<value_type> values(values_size, 0);
    std::generate(values.begin(), values.end(), [&]()
    {
        return value_distribution(random_engine);
    });

    return values;
}()};

template <auto transpose_matrix_fn, typename simd_t, size_t rows, size_t columns>
void transpose_matrix_benchmark(benchmark::State & state)
{
    constexpr size_t chunk_size = columns * rows;
    static_assert(simd_length<simd_t> == columns);

    std::vector<int, seqan3::aligned_allocator<int, sizeof(simd_t)>> values(_g_values.begin(), _g_values.end());
    simd_t sum{};

    simd_t * const matrix_begin = reinterpret_cast<simd_t *>(values.data());
    simd_t * const matrix_end = reinterpret_cast<simd_t *>(values.data() + ((values.size() / chunk_size) * chunk_size));

    // alignment is 128 bits / 8 bits/byte = 16 byte, 256 bits / 8 bits/byte = 32 byte
    assert((((size_t)matrix_begin) % sizeof(simd_t)) == 0);

    // ensure that each chunk has size of exactly e.g. 8x8 = 64
    assert(((reinterpret_cast<int32_t *>(matrix_end) - reinterpret_cast<int32_t *>(matrix_begin)) % (columns * rows)) == 0);

    for (auto _ : state)
    {
        for (simd_t * matrix_it = matrix_begin; matrix_it < matrix_end; matrix_it += rows)
        {
            std::span<simd_t, rows> matrix_span(matrix_it, rows);
            transpose_matrix_fn(matrix_span);

            for (size_t i = 0; i < rows; ++i)
                sum += matrix_span[i];
        }
    }
}

BENCHMARK_TEMPLATE(transpose_matrix_benchmark, transpose_matrix_32x4x4_omp, int32x4_t, 4, 4);
BENCHMARK_TEMPLATE(transpose_matrix_benchmark, transpose_matrix_32x4x4_sse4, int32x4_t, 4, 4);
BENCHMARK_TEMPLATE(transpose_matrix_benchmark, transpose_matrix_32x4x4_avx2, int32x4_t, 4, 4);
BENCHMARK_TEMPLATE(transpose_matrix_benchmark, transpose_matrix_32x4x4_avx2_gather, int32x4_t, 4, 4);

BENCHMARK_TEMPLATE(transpose_matrix_benchmark, transpose_matrix_32x8x8_omp, int32x8_t, 8, 8);
// BENCHMARK_TEMPLATE(transpose_matrix_benchmark, transpose_matrix_32x8x8_sse4, int32x8_t, 8, 8);
BENCHMARK_TEMPLATE(transpose_matrix_benchmark, transpose_matrix_32x8x8_avx2, int32x8_t, 8, 8);
BENCHMARK_TEMPLATE(transpose_matrix_benchmark, transpose_matrix_32x8x8_avx2_gather, int32x8_t, 8, 8);

BENCHMARK_TEMPLATE(transpose_matrix_benchmark, transpose_matrix_32x6x4_omp, int32x4_t, 6, 4);
BENCHMARK_TEMPLATE(transpose_matrix_benchmark, transpose_matrix_32x6x4_sse4, int32x4_t, 6, 4);
BENCHMARK_TEMPLATE(transpose_matrix_benchmark, transpose_matrix_32x6x4_avx2, int32x4_t, 6, 4);
BENCHMARK_TEMPLATE(transpose_matrix_benchmark, transpose_matrix_32x6x4_avx2_gather, int32x4_t, 6, 4);
