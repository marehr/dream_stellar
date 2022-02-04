#include <gtest/gtest.h>

#include <seqan3/std/span>

#include <stellar/prefilter/detail/valik_local_prefilter.hpp>
#include <stellar/prefilter/detail/valik_local_query_prefilter_index.hpp>

#include <stellar/test/fixture/512_simSeq_5e-2.hpp>

auto build_lqpindex()
{
    stellar::detail::valik_local_query_prefilter_index lqpindex
    {
        // ibf parameters
        ._ibf {
            seqan3::bin_count{5u /*arguments->bins*/},
            seqan3::bin_size{4096u / 1u/*arguments->bits / arguments->parts*/},
            seqan3::hash_function_count{2u/*arguments->hash*/}
        },

        // hash parameters
        ._kmer_size{seqan3::ungapped{4u/*arguments->kmer_size*/}},
        ._window_size{seqan3::window_size{6u/*arguments->window_size*/}}
    };

    return lqpindex;
}

TEST(valik_local_prefilter, index)
{
    auto lqpindex = build_lqpindex();
}

TEST(valik_local_prefilter, min_length_is_smaller_than_window_size)
{
    auto lqpindex = build_lqpindex();

    seqan::String<seqan::Dna5> sequence0{"TAGAC" "TAACCGCAGAACACG" "ACTCCTCTAC"};
    seqan::String<seqan::Dna5> sequence1{"TCACGCCAGTTTAATTTTAGGGCTCCAAGT"};
    seqan::String<seqan::Dna5> sequence2{"TACGCATATCTGG" "TAACCGCAGAACACG" "AA"};

    lqpindex.insert_sequence(seqan3::bin_index{0u}, sequence0);
    lqpindex.insert_sequence(seqan3::bin_index{1u}, sequence1);
    lqpindex.insert_sequence(seqan3::bin_index{2u}, sequence2);

    size_t const minLength = 5u;
    size_t const maxError = 0u;
    seqan::String<seqan::Dna5> raw_query{"TAACCGCAGAACACG"};

    for (size_t query_length = minLength; query_length < 15u; ++query_length)
    {
        seqan::String<seqan::Dna5> query = seqan::infix(raw_query, 0u, query_length);

        // Execute query on all databases, as the prefilter can't determine if it is a fit
        auto on_result = [expected_bin_id = 0](size_t const bin_id) mutable
        {
            EXPECT_EQ(expected_bin_id, bin_id);
            ++expected_bin_id;
        };

        testing::internal::CaptureStderr();
        stellar::detail::valik_local_prefilter<seqan::Dna5>(query, lqpindex, minLength, maxError, on_result);

        std::string expected_error =
            "WARNING: minLength (= 5) is smaller than the windowSize (= 6) of the minimiser. "
            "The prefilter can't be used and thus this query will be searched in all databases.\n";
        EXPECT_EQ(testing::internal::GetCapturedStderr(), expected_error);
    }
}

TEST(valik_local_prefilter, query_length_is_smaller_than_min_length)
{
    auto lqpindex = build_lqpindex();

    seqan::String<seqan::Dna5> sequence0{"TAGAC" "TAACCGCAGAACACG" "ACTCCTCTAC"};
    seqan::String<seqan::Dna5> sequence1{"TCACGCCAGTTTAATTTTAGGGCTCCAAGT"};
    seqan::String<seqan::Dna5> sequence2{"TACGCATATCTGG" "TAACCGCAGAACACG" "AA"};

    lqpindex.insert_sequence(seqan3::bin_index{0u}, sequence0);
    lqpindex.insert_sequence(seqan3::bin_index{1u}, sequence1);
    lqpindex.insert_sequence(seqan3::bin_index{2u}, sequence2);

    // query must be at least minLength long
    size_t const minLength = 15u;
    size_t const maxError = 0u;
    seqan::String<seqan::Dna5> raw_query{"TAACCGCAGAACACG"};

    for (size_t query_length = 0; query_length < 15u; ++query_length)
    {
        seqan::String<seqan::Dna5> query = seqan::infix(raw_query, 0u, query_length);

        // No results for queries that are shorter than min length
        auto on_result = [](size_t const bin_id)
        {
            FAIL() << "shouldn't be called, but was with bin_id: " << bin_id;
        };

        testing::internal::CaptureStderr();
        stellar::detail::valik_local_prefilter<seqan::Dna5>(query, lqpindex, minLength, maxError, on_result);

        std::string expected_error =
            "WARNING: sequence length (= " + std::to_string(query_length) + ") "
            "is smaller than the minLength (= 15). Skipping this sequence.\n";
        EXPECT_EQ(testing::internal::GetCapturedStderr(), expected_error);
    }
}

TEST(valik_local_prefilter, query_length_equals_min_length)
{
    auto lqpindex = build_lqpindex();

    seqan::String<seqan::Dna5> sequence0{"TAGAC" "TAACCGCAGAACACG" "ACTCCTCTAC"};
    seqan::String<seqan::Dna5> sequence1{"TCACGCCAGTTTAATTTTAGGGCTCCAAGT"};
    seqan::String<seqan::Dna5> sequence2{"TACGCATATCTGG" "TAACCGCAGAACACG" "AA"};

    lqpindex.insert_sequence(seqan3::bin_index{0u}, sequence0);
    lqpindex.insert_sequence(seqan3::bin_index{1u}, sequence1);
    lqpindex.insert_sequence(seqan3::bin_index{2u}, sequence2);

    // query must be at least minLength long
    size_t const minLength = 15u;
    size_t const maxError = 0u;
    seqan::String<seqan::Dna5> query0{seqan::infix(sequence0, 5u, 5u + minLength)};
    EXPECT_EQ(query0, "TAACCGCAGAACACG");

    std::vector<size_t> expected_bin_ids{0u, 2u};
    std::vector<size_t> bin_ids{};
    stellar::detail::valik_local_prefilter<seqan::Dna5>(query0, lqpindex, minLength, maxError, [&](size_t const bin_id)
    {
        bin_ids.push_back(bin_id);
    });

    EXPECT_EQ(bin_ids, expected_bin_ids);
}

TEST(valik_local_prefilter, query_length_is_longer_than_min_length)
{
    auto lqpindex = build_lqpindex();

    seqan::String<seqan::Dna5> sequence0{"TAGAC" "TAACCGCAGAACACG" "ACTCCTCTAC"};
    seqan::String<seqan::Dna5> sequence1{"TCACGCCAGTTTAATTTTAGGGCTCCAAGT"};
    seqan::String<seqan::Dna5> sequence2{"TACGCATATCTGG" "TAACCGCAGAACACG" "AA"};

    lqpindex.insert_sequence(seqan3::bin_index{0u}, sequence0);
    lqpindex.insert_sequence(seqan3::bin_index{1u}, sequence1);
    lqpindex.insert_sequence(seqan3::bin_index{2u}, sequence2);

    // query must be at least minLength long
    size_t const minLength = 15u;
    size_t const maxError = 0u;
    seqan::String<seqan::Dna5> query0{"TTTAATT" "TAACCGCAGAACACG" "CAAGT"};

    std::vector<size_t> expected_bin_ids{0u, 2u};
    std::vector<size_t> bin_ids{};
    stellar::detail::valik_local_prefilter<seqan::Dna5>(query0, lqpindex, minLength, maxError, [&](size_t const bin_id)
    {
        bin_ids.push_back(bin_id);
    });

    EXPECT_EQ(bin_ids, expected_bin_ids);
}

TEST(valik_local_prefilter, best_match_in_dataset_512_simSeq1_5eMinus2_fa)
{
    using TAlphabet = seqan::Dna5;
    stellar::test::alignment_fragment<TAlphabet> const & fragment
        = stellar::test::fixture::a512_simSeq_0dot05::forward::fragment0;

    stellar::detail::valik_local_query_prefilter_index lqpindex
    {
        // ibf parameters
        ._ibf {
            seqan3::bin_count{5u},
            seqan3::bin_size{524'288'000u / 16u}, /* = 32'768'000u*/
            // 32MB for 1MB sequence file LOL
            seqan3::hash_function_count{2u/*arguments->hash*/}
        },

        // hash parameters
        ._kmer_size{seqan3::ungapped{19u/*arguments->kmer_size*/}},
        ._window_size{seqan3::window_size{24u/*arguments->window_size*/}}
    };

    seqan::String<TAlphabet> sequence0{fragment.sequence1};
    lqpindex.insert_sequence(seqan3::bin_index{0u}, sequence0);

    // query must be at least minLength long
    size_t const minLength = 150u;
    double const errorRate = 0.05; // ⌈150u * 0.05⌉ = ⌈7.5⌉ = 8
    size_t const maxError = 2u;

    EXPECT_LE(maxError, std::ceil(minLength * errorRate));
    seqan::String<TAlphabet> query0{fragment.sequence2};

    std::vector<size_t> expected_bin_ids{0u};
    std::vector<size_t> bin_ids{};
    stellar::detail::valik_local_prefilter<seqan::Dna5>(query0, lqpindex, minLength, maxError, [&](size_t const bin_id)
    {
        bin_ids.push_back(bin_id);
    });

    EXPECT_EQ(bin_ids, expected_bin_ids);
}
