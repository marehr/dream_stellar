#include <gtest/gtest.h>

#include <stellar/prefilter/detail/valik_local_prefilter.hpp>
#include <stellar/prefilter/detail/valik_local_query_prefilter_index.hpp>

#include "build_lqpindex.hpp"
#include "minimiser_counter_state.hpp"
#include "minimiser_thresholder.hpp"

TEST(valik_local_prefilter_minimiser_counts, query15_minLength_15_maxError_0)
{
    auto lqpindex = build_lqpindex();

    seqan::String<seqan::Dna5> sequence0{"TAGAC" "TAACCGCAGAACACG" "ACTCCTCTAC"};
    seqan::String<seqan::Dna5> sequence1{"TCACGCCAGTTTAATTTTAGGGCTCCAAGT"};
    seqan::String<seqan::Dna5> sequence2{"TACGCATATCTGG" "TAACCGCAGAACACG" "AA"};

    lqpindex.insert_sequence(seqan3::bin_index{0u}, sequence0);
    lqpindex.insert_sequence(seqan3::bin_index{1u}, sequence1);
    lqpindex.insert_sequence(seqan3::bin_index{2u}, sequence2);

    size_t const minLength = 15u;
    size_t const maxError = 0u;

    seqan::String<seqan::Dna5> raw_query{"TAACCGCAGAACACG"};

    minimiser_counter_state counter{
        .lqpindex = lqpindex,
        .query_sequence = raw_query,
        .subquery_size = minLength
    };

    minimiser_thresholder thresholder{lqpindex, minLength, maxError};
    EXPECT_EQ(thresholder.min_number_of_minimisers(), 4u);
    EXPECT_EQ(thresholder.max_number_of_minimisers(), 10u);

    EXPECT_EQ(counter.subquery(), "TAACCGCAGAACACG");
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{32, 21, 29, 14, 52}));
    EXPECT_EQ(counter.minimiser_count(), 5u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 4u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{5, 0, 5, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{1, 0, 1, 0, 0}));
    EXPECT_FALSE(counter.next());
}

TEST(valik_local_prefilter_minimiser_counts, query27_minLength_15_maxError_0)
{
    auto lqpindex = build_lqpindex();

    seqan::String<seqan::Dna5> sequence0{"TAGAC" "TAACCGCAGAACACG" "ACTCCTCTAC"};
    seqan::String<seqan::Dna5> sequence1{"TCACGCCAGTTTAATTTTAGGGCTCCAAGT"};
    seqan::String<seqan::Dna5> sequence2{"TACGCATATCTGG" "TAACCGCAGAACACG" "AA"};

    lqpindex.insert_sequence(seqan3::bin_index{0u}, sequence0);
    lqpindex.insert_sequence(seqan3::bin_index{1u}, sequence1);
    lqpindex.insert_sequence(seqan3::bin_index{2u}, sequence2);

    size_t const minLength = 15u;
    size_t const maxError = 0u;

    seqan::String<seqan::Dna5> raw_query{"TTTAATT" "TAACCGCAGAACACG" "CAAGT"};

    minimiser_counter_state counter{
        .lqpindex = lqpindex,
        .query_sequence = raw_query,
        .subquery_size = minLength
    };

    minimiser_thresholder thresholder{lqpindex, minLength, maxError};
    EXPECT_EQ(thresholder.min_number_of_minimisers(), 4u);
    EXPECT_EQ(thresholder.max_number_of_minimisers(), 10u);

    EXPECT_EQ(counter.subquery(), "TTTAATTTAACCGCA" /*"GAACACGCAAGT"*/);
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{76, 79, 51, 32, 21}));
    EXPECT_EQ(counter.minimiser_count(), 5u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 4u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{3, 2, 2, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{0, 0, 0, 0, 0}));

    // same minimiser
    EXPECT_EQ(counter.subquery(), /*"T"*/ "TTAATTTAACCGCAG" /*"AACACGCAAGT"*/);
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{76, 79, 51, 32, 21}));
    EXPECT_EQ(counter.minimiser_count(), 5u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 4u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{3, 2, 2, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{0, 0, 0, 0, 0}));

    // new minimiser
    EXPECT_EQ(counter.subquery(), /*"TT"*/ "TAATTTAACCGCAGA" /*"ACACGCAAGT"*/);
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{76, 79, 51, 32, 21, 29}));
    EXPECT_EQ(counter.minimiser_count(), 6u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 5u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{4, 2, 3, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{0, 0, 0, 0, 0}));

    // new minimiser
    EXPECT_EQ(counter.subquery(), /*"TTT"*/ "AATTTAACCGCAGAA" /*"CACGCAAGT"*/);
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{79, 51, 32, 21, 29}));
    EXPECT_EQ(counter.minimiser_count(), 5u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 4u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{4, 1, 3, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{1, 0, 0, 0, 0}));

    // new minimiser
    EXPECT_EQ(counter.subquery(), /*"TTTA"*/ "ATTTAACCGCAGAAC" /*"ACGCAAGT"*/);
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{79, 51, 32, 21, 29, 14}));
    EXPECT_EQ(counter.minimiser_count(), 6u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 5u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{5, 1, 4, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{1, 0, 0, 0, 0}));

    // new minimiser
    EXPECT_EQ(counter.subquery(), /*"TTTAA"*/ "TTTAACCGCAGAACA" /*"CGCAAGT"*/);
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{51, 32, 21, 29, 14}));
    EXPECT_EQ(counter.minimiser_count(), 5u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 4u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{5, 0, 4, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{1, 0, 1, 0, 0}));

    // new minimiser
    EXPECT_EQ(counter.subquery(), /*"TTTAAT"*/ "TTAACCGCAGAACAC" /*"GCAAGT"*/);
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{32, 21, 29, 14}));
    EXPECT_EQ(counter.minimiser_count(), 4u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 3u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{4, 0, 4, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{1, 0, 1, 0, 0}));

    // new minimiser
    EXPECT_EQ(counter.subquery(), /*"TTTAATT"*/ "TAACCGCAGAACACG" /*"CAAGT"*/);
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{32, 21, 29, 14, 52}));
    EXPECT_EQ(counter.minimiser_count(), 5u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 4u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{5, 0, 5, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{1, 0, 1, 0, 0}));

    // new minimiser
    EXPECT_EQ(counter.subquery(), /*"TTTAATTT"*/ "AACCGCAGAACACGC" /*"AAGT"*/);
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{21, 29, 14, 52, 20}));
    EXPECT_EQ(counter.minimiser_count(), 5u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 4u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{4, 1, 5, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{1, 0, 1, 0, 0}));

    // same minimiser
    EXPECT_EQ(counter.subquery(), /*"TTTAATTTA"*/ "ACCGCAGAACACGCA" /*"AGT"*/);
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{21, 29, 14, 52, 20}));
    EXPECT_EQ(counter.minimiser_count(), 5u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 4u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{4, 1, 5, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{1, 0, 1, 0, 0}));

    // same minimiser
    EXPECT_EQ(counter.subquery(), /*"TTTAATTTAA"*/ "CCGCAGAACACGCAA" /*"GT"*/);
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{21, 29, 14, 52, 20}));
    EXPECT_EQ(counter.minimiser_count(), 5u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 4u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{4, 1, 5, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{1, 0, 1, 0, 0}));

    // new minimiser
    EXPECT_EQ(counter.subquery(), /*"TTTAATTTAAC"*/ "CGCAGAACACGCAAG" /*"T"*/);
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{29, 14, 52, 20, 31}));
    EXPECT_EQ(counter.minimiser_count(), 5u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 4u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{3, 1, 4, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{0, 0, 1, 0, 0}));

    // same minimiser
    EXPECT_EQ(counter.subquery(), /*"TTTAATTTAACC"*/ "GCAGAACACGCAAGT");
    EXPECT_TRUE(counter.next());
    EXPECT_EQ(counter.minimiser_hashes(), (std::vector<size_t>{29, 14, 52, 20, 31}));
    EXPECT_EQ(counter.minimiser_count(), 5u);
    EXPECT_EQ(thresholder.get(counter.minimiser_count()), 4u);
    EXPECT_EQ(counter.bin_counts(),
              (std::vector<uint16_t>{3, 1, 4, 0, 0}));
    EXPECT_EQ(thresholder.bin_thresholds_met(counter.bin_counts(), counter.minimiser_count()),
              (std::vector<bool>{0, 0, 1, 0, 0}));
}
