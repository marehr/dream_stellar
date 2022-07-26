#include <benchmark/benchmark.h>

#include <random>

#include <stellar/test/prefilter/local_query_prefilter/minimiser_window/seqan3_minimiser_window.hpp>
#include <stellar/test/prefilter/local_query_prefilter/minimiser_window/stellar_minimiser_window.hpp>

using value_type = int;
static constexpr size_t values_size = 100'000;

std::vector<value_type> values{[]()
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

template <typename minimiser_window_t>
void minimiser_window(benchmark::State & state)
{
    static constexpr size_t window_size = 10;

    value_type minimiser_sum{};
    for (auto _ : state)
    {
        minimiser_window_t minimiser_window{window_size};

        auto it = minimiser_window.initialize(values.begin(), values.end());
        minimiser_sum += minimiser_window.min();

        for (auto sentinel = values.end(); it != sentinel; ++it)
        {
            bool new_minimiser = minimiser_window.cyclic_push(*it);
            if (new_minimiser)
                minimiser_sum += minimiser_window.min();
        }
    }
}

BENCHMARK_TEMPLATE(minimiser_window, stellar::test::seqan3_minimiser_window<int>);
BENCHMARK_TEMPLATE(minimiser_window, stellar::test::stellar_minimiser_window<int, false>);
BENCHMARK_TEMPLATE(minimiser_window, stellar::test::stellar_minimiser_window<int, true>);
