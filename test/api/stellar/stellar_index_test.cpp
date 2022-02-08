#include <gtest/gtest.h>

#include <stellar/stellar_index.hpp>
#include <stellar/stellar_database_segment.hpp>
#include <stellar/stellar_query_segment.hpp>
#include <stellar/stellar_query_segment.tpp>

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
};

using ConstructFactories = ::testing::Types<StringSetOwnerFactory, StringSetDependentFactory, InfixSegmentVectorFactory>;
TYPED_TEST_SUITE(StellarIndexTest, ConstructFactories);

TYPED_TEST(StellarIndexTest, Construct)
{
    using TAlphabet = seqan::Dna5;

    using TSequence = seqan::String<TAlphabet>;
    using TInfixSequence = seqan::Segment<TSequence const, seqan::InfixSegment>;

    TSequence database{"TAGAC" "TAACCGCAGAACACG" "ACTCCTCTAC"};

    seqan::StringSet<TSequence> queries{};
    seqan::appendValue(queries, "TACGCATATCTGG" "TAACCGCAGAACACG" "AA");
    seqan::appendValue(queries, "TAGCCG" "TAACCGCAGAACACG" "GACC");
    seqan::appendValue(queries, "TAGCCAGTTTA" "GCAGAACAC" "CAAGA");

    size_t const kmerSize = 7u;
    stellar::StellarOptions options{};
    options.qGram = kmerSize;
    options.epsilon = 0.01;
    options.minLength = 9u;

    decltype(auto) dependentQueries = this->factory.dependentQueries(queries);

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
    EXPECT_EQ(std::addressof(this->factory.underlyingSequence(dependentQueries[0])), std::addressof(queries[0]));
    EXPECT_EQ(std::addressof(this->factory.underlyingSequence(dependentQueries[1])), std::addressof(queries[1]));
    EXPECT_EQ(std::addressof(this->factory.underlyingSequence(dependentQueries[2])), std::addressof(queries[2]));

    stellar::StellarIndex<TAlphabet> index{dependentQueries, options};
    index.construct();

    stellar::StellarSwiftPattern<TAlphabet> swiftPattern = index.createSwiftPattern();
    stellar::StellarSwiftFinder<TAlphabet> swiftFinder
    {
        seqan::infix(database, 0u, seqan::length(database)),
        options.minRepeatLength,
        options.maxRepeatPeriod
    };

    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    using TQuerySegment = stellar::StellarQuerySegment<TAlphabet>;
    using TAlignment = std::pair<TDatabaseSegment, TQuerySegment>;

    std::vector<TAlignment> results{};

    while (seqan::find(swiftFinder, swiftPattern, options.epsilon, options.minLength)) {
        TDatabaseSegment const databaseSegment = TDatabaseSegment::fromFinderMatch(seqan::infix(swiftFinder));
        TQuerySegment const querySegment = TQuerySegment::fromPatternMatch(swiftPattern);

        results.emplace_back(databaseSegment, querySegment);
    }

    std::vector<std::pair<TDatabaseSegment, TQuerySegment>> expected{};
    using TInterval = std::pair<size_t, size_t>;

    ASSERT_EQ(results.size(), 3u);
    EXPECT_EQ(results[0].first.asInfixSegment(), "TAACCGCAGAACACGA");
    EXPECT_EQ(results[0].first.interval(), TInterval(5u, 5u + 16u));
    EXPECT_EQ(results[0].second.asInfixSegment(), "ATATCTGGTAACCGCAGAACACGAA");
    EXPECT_EQ(results[0].second.interval(), TInterval(5u, 5u + 25u));

    EXPECT_EQ(results[1].first.asInfixSegment(), "TAACCGCAGAACACG");
    EXPECT_EQ(results[1].first.interval(), TInterval(5u, 5u + 15u));
    EXPECT_EQ(results[1].second.asInfixSegment(), "GTAACCGCAGAACACGGACC");
    EXPECT_EQ(results[1].second.interval(), TInterval(5u, 5u + 20u));

    EXPECT_EQ(results[2].first.asInfixSegment(), "GCAGAACAC");
    EXPECT_EQ(results[2].first.interval(), TInterval(10u, 10u + 9u));
    EXPECT_EQ(results[2].second.asInfixSegment(), "AGCAGAACACCAAGA");
    EXPECT_EQ(results[2].second.interval(), TInterval(10u, 10u + 15u));
}
