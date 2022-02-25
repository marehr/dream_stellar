#include <gtest/gtest.h>

#include <stellar/stellar.hpp>

namespace seqan {

template <typename TAlphabet>
struct Cargo<::stellar::StellarQGramIndex<TAlphabet>>
{
    typedef struct
    {
        double      abundanceCut;
    } Type;
};

} // namespace seqan

struct StringSetOwnerFactory
{
    template <typename TAlphabet>
    static decltype(auto)
    dependentQueries(seqan::StringSet<seqan::String<TAlphabet>, seqan::Owner<> > const & sequences)
    {
        return sequences; // pass as is
    }

    template <typename TAlphabet>
    static decltype(auto) underlyingSequence(seqan::String<TAlphabet> const & sequence)
    {
        return sequence; // pass as is
    }
};

struct StringSetDependentFactory
{
    template <typename TAlphabet>
    static seqan::StringSet<seqan::String<TAlphabet>, seqan::Dependent<>>
    dependentQueries(seqan::StringSet<seqan::String<TAlphabet>> const & sequences)
    {
        return sequences; // convert
    }

    template <typename TAlphabet>
    static decltype(auto) underlyingSequence(seqan::String<TAlphabet> const & sequence)
    {
        return sequence; // pass as is
    }
};

struct InfixSegmentVectorFactory
{
    template <typename TAlphabet>
    using TInfixSequence = seqan::Segment<seqan::String<TAlphabet> const, seqan::InfixSegment>;

    template <typename TAlphabet>
    static auto
    dependentQueries(seqan::StringSet<seqan::String<TAlphabet>> const & sequences)
    {
        std::vector<TInfixSequence<TAlphabet>> infixSequences{};
        for (auto const & sequence : sequences)
            infixSequences.emplace_back(sequence, 0u, seqan::length(sequence));

        return infixSequences; // convert
    }

    template <typename TAlphabet>
    static decltype(auto) underlyingSequence(TInfixSequence<TAlphabet> const & sequence)
    {
        return host(sequence); // pass as is
    }
};

template <typename TFactory>
struct StellarIndexTest : public ::testing::Test {
public:
    TFactory factory{};

    template <typename TSequenceSegment, typename TSequence>
    static void expect_segment(
        TSequenceSegment const & sequenceSegment,
        TSequence const & expectedUnderlyingSequence,
        size_t const expectedBeginPosition,
        size_t const expectedEndPosition,
        std::string_view const expectedSubsequence)
    {
        using TInterval = std::pair<size_t, size_t>;

        TSequenceSegment expectedSequenceSegment{
            expectedUnderlyingSequence,
            expectedBeginPosition,
            expectedEndPosition
        };

        TSequence const & actualSequence = sequenceSegment.underlyingSequence();

        // expect that the infix segment points to the right underlying sequence
        EXPECT_EQ(std::addressof(actualSequence), std::addressof(expectedUnderlyingSequence));

        // check whether underlying sequence has the same content (it should if prior test passes)
        EXPECT_EQ(actualSequence, expectedUnderlyingSequence);

        // check whether sequence corresponds to expected sub-sequence
        EXPECT_EQ(sequenceSegment.asInfixSegment(), std::string(expectedSubsequence));

        // check whether segment interval is the same as expected
        EXPECT_EQ(sequenceSegment.interval(), TInterval(expectedBeginPosition, expectedEndPosition));

        // last check if is the same segment
        EXPECT_EQ(sequenceSegment, expectedSequenceSegment);
    }
};

using ConstructFactories = ::testing::Types<StringSetOwnerFactory, StringSetDependentFactory, InfixSegmentVectorFactory>;
TYPED_TEST_SUITE(StellarIndexTest, ConstructFactories);

template <typename TAlphabet>
struct results
{
    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    using TQuerySegment = stellar::StellarQuerySegment<TAlphabet>;
    using TAlignment = std::pair<TDatabaseSegment, TQuerySegment>;

    // Note: verifySwiftHit's delta is the sum of pattern_delta + pattern_overlap
    std::vector<unsigned> pattern_deltas{};
    std::vector<unsigned> pattern_overlaps{};
    std::vector<unsigned> verify_deltas{};
    std::vector<TAlignment> alignments{};
    stellar::StellarComputeStatistics statistics{};
};

// Hijack SwiftHitVerifier for testing.
namespace stellar
{
struct TVerifyPassthroughSeed{};

template <>
struct SwiftHitVerifier<TVerifyPassthroughSeed>
{
    double const epsilon;
    int const minLength;
    double const xDrop;

    template <typename TAlphabet, typename TDelta, typename TOnAlignmentResultFn>
    void verify(StellarDatabaseSegment<TAlphabet> const & databaseSegment,
                StellarQuerySegment<TAlphabet> const & querySegment,
                TDelta const delta,
                TOnAlignmentResultFn && onAlignmentResult)
    {
        onAlignmentResult(databaseSegment, querySegment, delta);
    }
};

} // namespace stellar

template <typename TAlphabet, typename TypeParam>
results<TAlphabet> stellar_kernel_swift_seeds(
    StellarIndexTest<TypeParam> const & self,
    seqan::String<TAlphabet> const & database,
    seqan::StringSet<seqan::String<TAlphabet>> const & queries,
    stellar::StellarOptions const & options)
{
    using TSequence = seqan::String<TAlphabet>;
    using TInfixSequence = seqan::Segment<TSequence const, seqan::InfixSegment>;

    decltype(auto) dependentQueries = self.factory.dependentQueries(queries);

    if constexpr (std::is_same<TypeParam, StringSetOwnerFactory>::value)
    {
        EXPECT_TRUE((std::is_same<seqan::StringSet<TSequence> const &, decltype(dependentQueries)>::value));
    } else if constexpr (std::is_same<TypeParam, StringSetDependentFactory>::value)
    {
        EXPECT_TRUE((std::is_same<seqan::StringSet<TSequence, seqan::Dependent<>>, decltype(dependentQueries)>::value));
    } else if constexpr (std::is_same<TypeParam, InfixSegmentVectorFactory>::value)
    {
        EXPECT_TRUE((std::is_same<std::vector<TInfixSequence>, decltype(dependentQueries)>::value));
    } else
    {
        ADD_FAILURE() << "Unknown Factory";
    }

    // make sure that dependentQueries are reference to the actual queries (have same pointer address)
    EXPECT_EQ(std::addressof(self.factory.underlyingSequence(dependentQueries[0])), std::addressof(queries[0]));
    EXPECT_EQ(std::addressof(self.factory.underlyingSequence(dependentQueries[1])), std::addressof(queries[1]));
    EXPECT_EQ(std::addressof(self.factory.underlyingSequence(dependentQueries[2])), std::addressof(queries[2]));

    stellar::StellarIndex<TAlphabet> index{dependentQueries, options};
    index.construct();

    stellar::StellarSwiftPattern<TAlphabet> swiftPattern = index.createSwiftPattern();
    stellar::StellarSwiftFinder<TAlphabet> swiftFinder
    {
        seqan::infix(database, 0u, seqan::length(database)),
        options.minRepeatLength,
        options.maxRepeatPeriod
    };

    static_assert(std::is_same<decltype(swiftPattern.bucketParams[0].delta), unsigned>::value, "");
    static_assert(std::is_same<decltype(swiftPattern.bucketParams[0].overlap), unsigned>::value, "");

    using TResults = results<TAlphabet>;
    using TDatabaseSegment = typename TResults::TDatabaseSegment;
    using TQuerySegment = typename TResults::TQuerySegment;

    TResults results{};

    stellar::SwiftHitVerifier<stellar::TVerifyPassthroughSeed> verifier
    {
        .epsilon = options.epsilon,
        .minLength = options.minLength,
        .xDrop = options.xDrop
    };
    auto isPatternDisabled = [](...){ return false; };

    results.statistics = stellar::_stellarKernel(
        swiftFinder,
        swiftPattern,
        verifier,
        isPatternDisabled,
        [&](TDatabaseSegment const & databaseSegment, TQuerySegment const & querySegment, unsigned const delta)
    {
        // should always be the same number
        results.pattern_deltas.emplace_back(swiftPattern.bucketParams[0].delta);
        // should always be the same number
        results.pattern_overlaps.emplace_back(swiftPattern.bucketParams[0].overlap);
        // should always be the same number = pattern_deltas + pattern_overlaps
        results.verify_deltas.emplace_back(delta);
        // different for each invocation
        results.alignments.emplace_back(databaseSegment, querySegment);
    });

    // sort alignments by query and not database
    std::sort(results.alignments.begin(), results.alignments.end(), [&](auto const & v1, auto const & v2)
    {
        return v1.second < v2.second;
    });

    return results;
}

TYPED_TEST(StellarIndexTest, validSeedWithExactMatches)
{
    using TAlphabet = seqan::Dna5;
    using TSequence = seqan::String<TAlphabet>;
    TSequence database{"CAACGGACTGCTGTCTAGAC" "TAACCGCAGAACACG" "ACTCCTCTACCTTACCGCGT"};

    seqan::StringSet<TSequence> queries{};
    // original
    seqan::appendValue(queries, "CTCGAGGGTTTACGCATATCTGG" "TAACCG" "C"/**/ "AGAACACG" "A" "AGAGCCTGAGA");
    // 1 error, insert in query sequence
    seqan::appendValue(queries, "CTCGAGGGTTTACGCATATCTGG" "TAACCG" "C" "A" "AGAACACG" "A" "AGAGCCTGAGA");
    // 1 error, delete in query sequence
    seqan::appendValue(queries, "CTCGAGGGTTTACGCATATCTGG" "TAACCG"/*C   */ "AGAACACG" "A" "AGAGCCTGAGA");
    // 1 error, substitution in query sequence
    seqan::appendValue(queries, "CTCGAGGGTTTACGCATATCTGG" "TAACCG" "G"/**/ "AGAACACG" "A" "AGAGCCTGAGA");
    // another sequence
    seqan::appendValue(queries, "CCGACTACCCACTTACTTATTAGCCG" "TAACCGCAGAACACG" "GACCAATCAGGCCC");
    // another smaller sequence
    seqan::appendValue(queries, "TAGCCAGTTTA" "GCAGAACAC" "CAAGA");

    size_t const kmerSize = 7u;
    stellar::StellarOptions options{};
    options.qGram = kmerSize;
    options.epsilon = 0.01;
    options.minLength = 9u;

    using TResults = results<TAlphabet>;
    TResults results = stellar_kernel_swift_seeds(*this, database, queries, options);

    EXPECT_EQ(results.statistics.numSwiftHits, 6u);
    EXPECT_EQ(results.statistics.totalLength, 16u + 9u + 9u + 9u + 15u + 9u);
    EXPECT_EQ(results.statistics.maxLength, 16u);
    EXPECT_EQ(results.verify_deltas, (std::vector<unsigned>{16u, 16u, 16u, 16u, 16u, 16u}));
    EXPECT_EQ(results.pattern_deltas, (std::vector<unsigned>{16u, 16u, 16u, 16u, 16u, 16u}));
    EXPECT_EQ(results.pattern_overlaps, (std::vector<unsigned>{0u, 0u, 0u, 0u, 0u, 0u}));

    ASSERT_EQ(results.alignments.size(), 6u);

    // original
    this->expect_segment(results.alignments[0].first, database, 20u, 20u + 16u,
                         /*  */"TAACCGCAGAACACGA");
    this->expect_segment(results.alignments[0].second, queries[0], 20u, 20u + 3u + 16u + 11u,
                         "TGG" "TAACCGCAGAACACGA" "AGAGCCTGAGA");

    // 1 error, insert in query sequence
    this->expect_segment(results.alignments[1].first, database, 27u, 27u + 9u,
                         /*   */"AGAACACGA");
    this->expect_segment(results.alignments[1].second, queries[1], 27u, 27u + 4u + 9u + 11u,
                         "CGCA" "AGAACACGA" "AGAGCCTGAGA");

    // 1 error, delete in query sequence
    this->expect_segment(results.alignments[2].first, database, 27u, 27u + 9u,
                         /* */"AGAACACGA");
    this->expect_segment(results.alignments[2].second, queries[2], 27u, 27u + 2u + 9u + 11u,
                         "CG" "AGAACACGA" "AGAGCCTGAGA");

    // 1 error, substitution in query sequence
    this->expect_segment(results.alignments[3].first, database, 27u, 27u + 9u,
                         /*  */"AGAACACGA");
    this->expect_segment(results.alignments[3].second, queries[3], 27u, 27u + 3u + 9u + 11u,
                         "CGG" "AGAACACGA" "AGAGCCTGAGA");

    // another sequence
    this->expect_segment(results.alignments[4].first, database, 20u, 20u + 15u,
                         /*     */"TAACCGCAGAACACG");
    this->expect_segment(results.alignments[4].second, queries[4], 20u, 20u + 6u + 15u + 10u,
                         "TAGCCG" "TAACCGCAGAACACG" "GACCAATCAG");

    // another smaller sequence
    this->expect_segment(results.alignments[5].first, database, 25u, 25u + 9u,
                         /* */"GCAGAACAC");
    this->expect_segment(results.alignments[5].second, queries[5], 9u, 9u + 2u + 9u + 5u,
                         "TA" "GCAGAACAC" "CAAGA");
}

TYPED_TEST(StellarIndexTest, validSeedWithOneErrorMatches)
{
    using TAlphabet = seqan::Dna5;
    using TSequence = seqan::String<TAlphabet>;
    TSequence database{"CAACGGACTGCTGTCTAGAC" "TAACCGCAGAACACG" "ACTCCTCTACCTTACCGCGT"};

    seqan::StringSet<TSequence> queries{};
    // original
    seqan::appendValue(queries, "CTCGAGGGTTTACGCATATCTGG" "TAACCG" "C"/**/ "AGAACACG" "A" "AGAGCCTGAGA");
    // 1 error, insert in query sequence
    seqan::appendValue(queries, "CTCGAGGGTTTACGCATATCTGG" "TAACCG" "C" "A" "AGAACACG" "A" "AGAGCCTGAGA");
    // 1 error, delete in query sequence
    seqan::appendValue(queries, "CTCGAGGGTTTACGCATATCTGG" "TAACCG"/*C   */ "AGAACACG" "A" "AGAGCCTGAGA");
    // 1 error, substitution in query sequence
    seqan::appendValue(queries, "CTCGAGGGTTTACGCATATCTGG" "TAACCG" "G"/**/ "AGAACACG" "A" "AGAGCCTGAGA");
    // another sequence
    seqan::appendValue(queries, "CCGACTACCCACTTACTTATTAGCCG" "TAACCGCAGAACACG" "GACCAATCAGGCCC");
    // another smaller sequence
    seqan::appendValue(queries, "TAGCCAGTTTA" "GCAGAACAC" "CAAGA");

    size_t const kmerSize = 7u;
    stellar::StellarOptions options{};
    options.qGram = kmerSize;
    options.epsilon = 0.125; // 9 length * 0.125 error-rate = 1 error
    options.minLength = 9u;

    using TResults = results<TAlphabet>;
    TResults results = stellar_kernel_swift_seeds(*this, database, queries, options);

    EXPECT_EQ(results.statistics.numSwiftHits, 6u);
    EXPECT_EQ(results.statistics.totalLength, 16u + 16u + 9u + 9u + 15u + 9u);
    EXPECT_EQ(results.statistics.maxLength, 16u);
    EXPECT_EQ(results.verify_deltas, (std::vector<unsigned>{22u, 22u, 22u, 22u, 22u, 22u}));
    EXPECT_EQ(results.pattern_deltas, (std::vector<unsigned>{16u, 16u, 16u, 16u, 16u, 16u}));
    EXPECT_EQ(results.pattern_overlaps, (std::vector<unsigned>{6u, 6u, 6u, 6u, 6u, 6u}));

    ASSERT_EQ(results.alignments.size(), 6u);
    // original (same result as with 0 errors, but 6 chars earlier)
    this->expect_segment(results.alignments[0].first, database, 20u, 20u + 16u,
                         /*           */"TAACCGCAGAACACGA");
    this->expect_segment(results.alignments[0].second, queries[0], 14u, 14u + 6u + 3u + 16u + 11u,
                         "CATATC" "TGG" "TAACCGCAGAACACGA" "AGAGCCTGAGA");

    // 1 error, insert in query sequence (complete new result compared with 0 errors)
    this->expect_segment(results.alignments[1].first, database, 20u, 20u + 16u,
                         /*        */"TAACCGC"/* */"AGAACACGA");
    this->expect_segment(results.alignments[1].second, queries[1], 14u, 14u + 9u + 17u + 11u,
                         "CATATCTGG" "TAACCGC" "A" "AGAACACGA" "AGAGCCTGAGA");

    // 1 error, delete in query sequence (same result as with 0 errors, but 6 chars earlier)
    this->expect_segment(results.alignments[2].first, database, 27u, 27u + 9u,
                         /*          */"AGAACACGA");
    this->expect_segment(results.alignments[2].second, queries[2], 21u, 21u + 6u + 2u + 9u + 11u,
                         "GGTAAC" "CG" "AGAACACGA" "AGAGCCTGAGA");

    // 1 error, substitution in query sequence (same result as with 0 errors, but 6 chars earlier)
    this->expect_segment(results.alignments[3].first, database, 27u, 27u + 9u,
                         /*          */"AGAACACGA");
    this->expect_segment(results.alignments[3].second, queries[3], 21u, 21u + 6u + 3u + 9u + 11u,
                        "GGTAAC" "CGG" "AGAACACGA" "AGAGCCTGAGA");

    // another sequence (same result as with 0 errors, but 6 chars earlier)
    this->expect_segment(results.alignments[4].first, database, 20u, 20u + 15u,
                         /*              */"TAACCGCAGAACACG");
    this->expect_segment(results.alignments[4].second, queries[4], 14u, 14u + 6u + 6u + 15u + 10u,
                         "ACTTAT" "TAGCCG" "TAACCGCAGAACACG" "GACCAATCAG");

    // another smaller sequence (same result as with 0 errors, but 6 chars earlier)
    this->expect_segment(results.alignments[5].first, database, 25u, 25u + 9u,
                         /*          */"GCAGAACAC");
    this->expect_segment(results.alignments[5].second, queries[5], 3u, 3u + 6u + 2u + 9u + 5u,
                         "CCAGTT" "TA" "GCAGAACAC" "CAAGA");
}
