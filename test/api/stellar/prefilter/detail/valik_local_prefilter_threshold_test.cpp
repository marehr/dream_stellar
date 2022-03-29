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
    raptor::threshold::threshold_parameters threshold_parameters = search_arguments.make_threshold_parameters();

    // kmers = 15 - 4 + 1 = 12
    std::vector<size_t> threshold_probabilistic_data = raptor::threshold::precompute_threshold(threshold_parameters);
    std::vector<size_t> threshold_probabilistic_data_expected{4, 5, 6, 7, 8, 9, 10};

    std::vector<size_t> threshold_probabilistic_correction = raptor::threshold::precompute_correction(threshold_parameters);
    std::vector<size_t> threshold_probabilistic_correction_expected{0, 0, 0, 0, 0, 0, 0};

    EXPECT_EQ(threshold_probabilistic_data, threshold_probabilistic_data_expected);
    EXPECT_EQ(threshold_probabilistic_correction, threshold_probabilistic_correction_expected);
}

TEST(valik_local_prefilter_threshold, minLength_150_maxError_2)
{
    auto lqpindex = build_lqpindex2();

    size_t const minLength = 150u;
    size_t const maxError = 2u;

    auto search_arguments = stellar::detail::valik_local_prefilter_search_arguments(lqpindex, minLength, maxError);
    raptor::threshold::threshold_parameters threshold_parameters = search_arguments.make_threshold_parameters();

    // kmers = 150 - 19 + 1 = 132
    std::vector<size_t> threshold_probabilistic_data = raptor::threshold::precompute_threshold(threshold_parameters);
    std::vector<size_t> threshold_probabilistic_data_expected{
#if 0
        // old threshold
          9,  10,  11,  11,  12,  12,  13,  14,  14,  15,  16,  16,  17,  18,  18,  19,  20,  20,  21,  22,  22,  23,
         24,  24,  25,  26,  27,  27,  28,  29,  29,  30,  31,  31,  32,  33,  34,  34,  35,  36,  36,  37,  38,  39,
         39,  40,  41,  41,  42,  43,  44,  44,  45,  46,  47,  47,  48,  49,  50,  50,  51,  52,  53,  53,  54,  55,
         56,  56,  57,  58,  59,  59,  60,  61,  62,  63,  63,  64,  65,  66,  67,  67,  68,  69,  70,  71,  71,  72,
         73,  74,  75,  76,  77,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89
#endif
          7,   8,   9,   9,  10,  11,  12,  12,  13,  14,  15,  15,  16,  17,  17,  18,  19,  20,  21,  21,  22,  23,
         24,  25,  26,  27,  28,  29,  30,  31,  31,  32,  33,  34,  35,  36,  37,  38,  39,  39,  40,  41,  42,  43,
         44,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,
         88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
        110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127
    };

    std::vector<size_t> threshold_probabilistic_correction = raptor::threshold::precompute_correction(threshold_parameters);
    std::vector<size_t> threshold_probabilistic_correction_expected(106, 0); // {0}x106

    EXPECT_EQ(threshold_probabilistic_data, threshold_probabilistic_data_expected);
    EXPECT_EQ(threshold_probabilistic_correction, threshold_probabilistic_correction_expected);
}
