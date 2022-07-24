
#include <gtest/gtest.h>

#include <stellar/test/prefilter/local_query_prefilter/minimiser_window/seqan3_minimiser_window.hpp>
#include <stellar/test/prefilter/local_query_prefilter/minimiser_window/stellar_minimiser_window.hpp>

#include <numeric>
#include <random>

template <typename TMinimiserWindow>
struct minimiser_window_test : public ::testing::Test
{};

using AgentSplitter = ::testing::Types<
    stellar::test::seqan3_minimiser_window<int>
    ,stellar::test::stellar_minimiser_window<int>
>;
TYPED_TEST_SUITE(minimiser_window_test, AgentSplitter);

TYPED_TEST(minimiser_window_test, increasing_sequence)
{
    std::vector<int> values{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    TypeParam minimiser_window{4};

    auto it = minimiser_window.initialize(values.begin(), values.end());

    // [0, 1, 2, 3], 4, 5, 6, 7, 8, 9, 10
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 3);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 0, [1, 2, 3, 4], 5, 6, 7, 8, 9, 10
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 4);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 1);
    ++it;

    // 0, 1, [2, 3, 4, 5], 6, 7, 8, 9, 10
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 5);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 2);
    ++it;

    // 0, 1, 2, [3, 4, 5, 6], 7, 8, 9, 10
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 6);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 3);
    ++it;

    // 0, 1, 2, 3, [4, 5, 6, 7], 8, 9, 10
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 7);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 4);
    ++it;

    // 0, 1, 2, 3, 4, [5, 6, 7, 8], 9, 10
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 8);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 5);
    ++it;

    // 0, 1, 2, 3, 4, 5, [6, 7, 8, 9], 10
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 9);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 6);
    ++it;

    // 0, 1, 2, 3, 4, 5, 6, [7, 8, 9, 10]
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 10);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 7);
    ++it;

    EXPECT_EQ(it, values.end());
}

TYPED_TEST(minimiser_window_test, decreasing_sequence)
{
    std::vector<int> values{10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

    TypeParam minimiser_window{4};

    auto it = minimiser_window.initialize(values.begin(), values.end());

    // [10, 9, 8, 7], 6, 5, 4, 3, 2, 1, 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 3);
    EXPECT_EQ(minimiser_window.min(), 7);
    ++it;

    // 10, [9, 8, 7, 6], 5, 4, 3, 2, 1, 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 4);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 6);
    ++it;

    // 10, 9, [8, 7, 6, 5], 4, 3, 2, 1, 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 5);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 5);
    ++it;

    // 10, 9, 8, [7, 6, 5, 4], 3, 2, 1, 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 6);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 4);
    ++it;

    // 10, 9, 8, 7, [6, 5, 4, 3], 2, 1, 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 7);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 3);
    ++it;

    // 10, 9, 8, 7, 6, [5, 4, 3, 2], 1, 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 8);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 2);
    ++it;

    // 10, 9, 8, 7, 6, 5, [4, 3, 2, 1], 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 9);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 1);
    ++it;

    // 10, 9, 8, 7, 6, 5, 4, [3, 2, 1, 0]
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 10);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    EXPECT_EQ(it, values.end());
}

TYPED_TEST(minimiser_window_test, random_sequence)
{
    std::vector<int> values{2, 4, 12, 0, 15, 13, 10, 9, 14, 8, 3, 1, 5, 11, 7, 6};

    TypeParam minimiser_window{4};

    auto it = minimiser_window.initialize(values.begin(), values.end());

    // [2, 4, 12, 0], 15, 13, 10, 9, 14, 8, 3, 1, 5, 11, 7, 6
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 3);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 2, [4, 12, 0, 15], 13, 10, 9, 14, 8, 3, 1, 5, 11, 7, 6
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 4);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 2, 4, [12, 0, 15, 13], 10, 9, 14, 8, 3, 1, 5, 11, 7, 6
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 5);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 2, 4, 12, [0, 15, 13, 10], 9, 14, 8, 3, 1, 5, 11, 7, 6
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 6);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 2, 4, 12, 0, [15, 13, 10, 9], 14, 8, 3, 1, 5, 11, 7, 6
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 7);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 9);
    ++it;

    // 2, 4, 12, 0, 15, [13, 10, 9, 14], 8, 3, 1, 5, 11, 7, 6
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 8);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 9);
    ++it;

    // 2, 4, 12, 0, 15, 13, [10, 9, 14, 8], 3, 1, 5, 11, 7, 6
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 9);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 8);
    ++it;

    // 2, 4, 12, 0, 15, 13, 10, [9, 14, 8, 3], 1, 5, 11, 7, 6
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 10);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 3);
    ++it;

    // 2, 4, 12, 0, 15, 13, 10, 9, [14, 8, 3, 1], 5, 11, 7, 6
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 11);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 1);
    ++it;

    // 2, 4, 12, 0, 15, 13, 10, 9, 14, [8, 3, 1, 5], 11, 7, 6
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 12);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 1);
    ++it;

    // 2, 4, 12, 0, 15, 13, 10, 9, 14, 8, [3, 1, 5, 11], 7, 6
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 13);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 1);
    ++it;

    // 2, 4, 12, 0, 15, 13, 10, 9, 14, 8, 3, [1, 5, 11, 7], 6
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 14);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 1);
    ++it;

    // 2, 4, 12, 0, 15, 13, 10, 9, 14, 8, 3, 1, [5, 11, 7, 6]
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 15);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 5);
    ++it;

    EXPECT_EQ(it, values.end());
}

TYPED_TEST(minimiser_window_test, same_sequence)
{
    std::vector<int> values{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    TypeParam minimiser_window{4};

    auto it = minimiser_window.initialize(values.begin(), values.end());

    // [0, 0, 0, *0], 0, 0, 0, 0, 0, 0, 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 3);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 0, [0, 0, *0, 0], 0, 0, 0, 0, 0, 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 4);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 0, 0, [0, *0, 0, 0], 0, 0, 0, 0, 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 5);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 0, 0, 0, [*0, 0, 0, 0], 0, 0, 0, 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 6);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 0, 0, 0, 0, [0, 0, 0, *0], 0, 0, 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 7);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 0, 0, 0, 0, 0, [0, 0, *0, 0], 0, 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 8);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 0, 0, 0, 0, 0, 0, [0, *0, 0, 0], 0
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 9);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 0, 0, 0, 0, 0, 0, 0, [*0, 0, 0, 0]
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 10);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    EXPECT_EQ(it, values.end());
}

TYPED_TEST(minimiser_window_test, regression01)
{
    std::vector<int> values{14, 8, 19, 3, 3, 9, 5, 20};

    TypeParam minimiser_window{3};

    auto it = minimiser_window.initialize(values.begin(), values.end());

    // [14, *8, 19], 3, 3, 9, 5, 20
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 2);
    EXPECT_EQ(minimiser_window.min(), 8);
    ++it;

    // 14, [8, 19, *3], 3, 9, 5, 20
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 3);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 3);
    ++it;

    // 14, 8, [19, *3, 3], 9, 5, 20
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 4);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 3);
    ++it;

    // 14, 8, 19, [*3, 3, 9], 5, 20
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 5);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 3);
    ++it;

    // 14, 8, 19, 3, [*3, 9, 5], 20
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 6);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 3);
    ++it;

    // 14, 8, 19, 3, 3, [9, *5, 20]
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 7);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 5);
    ++it;

    EXPECT_EQ(it, values.end());
}

TYPED_TEST(minimiser_window_test, regression02)
{
    std::vector<int> values{0, 1, 4, 0, 8, 18, 8, 16, 14, 9};

    TypeParam minimiser_window{3};

    auto it = minimiser_window.initialize(values.begin(), values.end());

    // [*0, 1, 4], 0, 8, 18, 8, 16, 14, 9
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 2);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 0, [1, 4, *0], 8, 18, 8, 16, 14, 9
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 3);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 0, 1, [4, *0, 8], 18, 8, 16, 14, 9
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 4);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 0, 1, 4, [*0, 8, 18], 8, 16, 14, 9
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 5);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 0);
    ++it;

    // 0, 1, 4, 0, [8, 18, *8], 16, 14, 9
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 6);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 8);
    ++it;

    // 0, 1, 4, 0, 8, [18, *8, 16], 14, 9
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 7);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 8);
    ++it;

    // 0, 1, 4, 0, 8, 18, [*8, 16, 14], 9
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 8);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), false);
    EXPECT_EQ(minimiser_window.min(), 8);
    ++it;

    // 0, 1, 4, 0, 8, 18, 8, [16, 14, *9]
    EXPECT_NE(it, values.end());
    EXPECT_EQ(it - values.begin(), 9);
    EXPECT_EQ(minimiser_window.cyclic_push(*it), true);
    EXPECT_EQ(minimiser_window.min(), 9);
    ++it;

    EXPECT_EQ(it, values.end());
}

TYPED_TEST(minimiser_window_test, random_reference_impl)
{
    std::random_device random_device{};
    std::size_t random_seed = 0;//random_device();
    std::mt19937_64 random_engine(random_seed);
    std::uniform_int_distribution<> distribution(0, 20);

    size_t window_size = std::uniform_int_distribution<>{3, 8}(random_engine);
    size_t values_size = std::uniform_int_distribution<>{20, 100}(random_engine);
    window_size = std::min(window_size, values_size);

    std::vector<int> values(values_size, 0);

    std::generate(values.begin(), values.end(), [&]()
    {
        return distribution(random_engine);
    });

    TypeParam minimiser_window{window_size};

    auto it = minimiser_window.initialize(values.begin(), values.end());

    auto value_window_it = values.begin();
    auto value_window_sentinel = value_window_it + window_size;

    auto minimiser = [](auto it, auto sen)
    {
        return std::min_element(it, sen, std::less_equal<int>{});
    };

    auto minimiser_it = minimiser(value_window_it, value_window_sentinel);
    EXPECT_EQ(minimiser_window.min(), *minimiser_it) << "seed: " << random_seed;

    auto diagnostics = [&]() -> std::string
    {
        std::string str;

        str += "seed: " + std::to_string(random_seed) + ", ";

        auto it = values.begin();

        for (; it != value_window_it; ++it)
            str += std::to_string(*it) + ", ";

        str += "[";
        for (; it != value_window_sentinel; ++it)
            str += (it == minimiser_it ? "!" : "") + std::to_string(*it) + ", ";
        str += "]";

        for (; it != values.end(); ++it)
            str += std::to_string(*it) + ", ";

        return str;
    };

    for (; value_window_sentinel != values.end();)
    {
        ++it;
        assert(it != values.end());
        assert(it == value_window_sentinel);
        ++value_window_it;
        ++value_window_sentinel;

        auto new_minimiser_it = minimiser(value_window_it, value_window_sentinel);
        bool minimiser_left_window = minimiser_it < value_window_it;
        bool smaller_overall_minimiser = *new_minimiser_it < *minimiser_it;
        if (minimiser_left_window || smaller_overall_minimiser)
            minimiser_it = new_minimiser_it;

        bool new_minimiser = minimiser_window.cyclic_push(*it);
        EXPECT_EQ(new_minimiser, minimiser_left_window || smaller_overall_minimiser) << diagnostics();

        EXPECT_EQ(minimiser_window.min(), *minimiser_it) << diagnostics();

    }
}
