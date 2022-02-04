#include <gtest/gtest.h>

#include <stellar/prefilter/detail/valik_local_prefilter_search_arguments.hpp>
#include <stellar/prefilter/detail/valik_local_prefilter_threshold.hpp>

#include "build_lqpindex.hpp"

TEST(valik_local_prefilter_threshold, minLength_15_maxError_0)
{
    auto lqpindex = build_lqpindex();

    size_t const minLength = 15u;
    size_t const maxError = 0u;

    auto search_arguments = stellar::detail::valik_local_prefilter_search_arguments(lqpindex, minLength, maxError);
    auto threshold = stellar::detail::valik_local_prefilter_make_threshold(search_arguments);

    std::vector<size_t> threshold_probabilistic_data{};
    std::vector<size_t> threshold_probabilistic_data_expected{4, 5, 6, 7, 8, 9, 10};

    std::vector<int> threshold_probabilistic_correction{};
    std::vector<int> threshold_probabilistic_correction_expected{-3, -3, -3, -3, -3, -3, -3};

    for (size_t i = 0; i < threshold.precomp_thresholds.size(); ++i)
    {
        size_t threshold_base = threshold.precomp_thresholds[i] + 3;
        int threshold_correction = (int)threshold.precomp_thresholds[i] - (int)threshold_base;
        threshold_probabilistic_data.push_back(threshold_base);
        threshold_probabilistic_correction.push_back(threshold_correction);
    }

    EXPECT_EQ(threshold_probabilistic_data, threshold_probabilistic_data_expected);
    EXPECT_EQ(threshold_probabilistic_correction, threshold_probabilistic_correction_expected);
}

TEST(valik_local_prefilter_threshold, minLength_150_maxError_2)
{
    auto lqpindex = build_lqpindex2();

    size_t const minLength = 150u;
    size_t const maxError = 2u;

    auto search_arguments = stellar::detail::valik_local_prefilter_search_arguments(lqpindex, minLength, maxError);
    auto threshold = stellar::detail::valik_local_prefilter_make_threshold(search_arguments);

    std::vector<size_t> threshold_probabilistic_data{};
    std::vector<size_t> threshold_probabilistic_data_expected{
        9, 10, 11, 11, 12, 12, 13, 14, 14, 15, 16, 16, 17, 18, 18, 19, 20, 20, 21, 22, 22, 23, 24, 24, 25, 26, 27, 27,
        28, 29, 29, 30, 31, 31, 32, 33, 34, 34, 35, 36, 36, 37, 38, 39, 39, 40, 41, 41, 42, 43, 44, 44, 45, 46, 47, 47,
        48, 49, 50, 50, 51, 52, 53, 53, 54, 55, 56, 56, 57, 58, 59, 59, 60, 61, 62, 63, 63, 64, 65, 66, 67, 67, 68, 69,
        70, 71, 71, 72, 73, 74, 75, 76, 77, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89
    };

    std::vector<int> threshold_probabilistic_correction{};
    std::vector<int> threshold_probabilistic_correction_expected(106,  -3);  // {-3}x106

    for (size_t i = 0; i < threshold.precomp_thresholds.size(); ++i)
    {
        size_t threshold_base = threshold.precomp_thresholds[i] + 3;
        int threshold_correction = (int)threshold.precomp_thresholds[i] - (int)threshold_base;
        threshold_probabilistic_data.push_back(threshold_base);
        threshold_probabilistic_correction.push_back(threshold_correction);
    }

    EXPECT_EQ(threshold_probabilistic_data, threshold_probabilistic_data_expected);
    EXPECT_EQ(threshold_probabilistic_correction, threshold_probabilistic_correction_expected);
}
