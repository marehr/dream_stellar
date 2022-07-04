#include <gtest/gtest.h>

#include <stellar/stellar_types.hpp>
#include <stellar/prefilter/no_query_prefilter.hpp>
#include <stellar/prefilter/whole_database_agent_splitter.hpp>
#include <stellar/prefilter/nsegment_database_agent_splitter.hpp>

template <typename TAgentSplitter>
struct NoQueryPrefilterTest : public ::testing::Test
{};

using AgentSplitter = ::testing::Types<stellar::WholeDatabaseAgentSplitter, stellar::NSegmentDatabaseAgentSplitter>;
TYPED_TEST_SUITE(NoQueryPrefilterTest, AgentSplitter);

TYPED_TEST(NoQueryPrefilterTest, prefilter)
{
    using TAlphabet = seqan::Dna5;
    using TAgentSplitter = TypeParam;

    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGGGGGGGGG");

    seqan::StringSet<seqan::String<TAlphabet> > queries{};
    appendValue(queries, "AAAA");
    appendValue(queries, "CCCC");
    appendValue(queries, "GGGG");
    appendValue(queries, "TTTT");

    using TPrefilter = stellar::NoQueryPrefilter<TAlphabet, TAgentSplitter>;
    using TPrefilterAgent = stellar::prefilter_agent<TAlphabet>;

    stellar::StellarOptions options{};
    options.minLength = 2u;
    options.threadCount = 3u;
    options.epsilon = {0, 1}; // 0.0

    TAgentSplitter const splitter = []()
    {
        if constexpr (std::is_same_v<TAgentSplitter, stellar::NSegmentDatabaseAgentSplitter>)
        {
            return stellar::NSegmentDatabaseAgentSplitter{4u};
        } else
        {
            return stellar::WholeDatabaseAgentSplitter{};
        }
    }();
    TPrefilter prefilter{options, databases, queries, splitter};

    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    std::vector<TDatabaseSegment> seenDatabases{};
    std::vector<TDatabaseSegment> expectedSeenDatabases = splitter.split(databases, options.minLength);

    // implicit conversion to base class
    for (TPrefilterAgent & prefilterAgent: prefilter.agents(options.threadCount, options))
    {
        using TDatabaseSegments = typename TPrefilter::TDatabaseSegments;
        using TQuerySwiftFilter = typename TPrefilter::TQueryFilter;

        prefilterAgent.prefilter([&](TDatabaseSegments databaseSegments, TQuerySwiftFilter & queryFilter)
        {
            // collect all seen databases
            for (TDatabaseSegment const & databaseSegment : databaseSegments)
                seenDatabases.push_back(databaseSegment);

            // This ensures that the queryFilter is been constructed with ALL queries
            auto && querySegments = stellar::StellarIndex<TAlphabet>::sequencesFromPattern(queryFilter);
            static_assert(std::is_same_v<decltype(querySegments), stellar::StellarQGramStringSet<TAlphabet> const &>);
            EXPECT_EQ(queries, querySegments);
        });
    }

    // Note: every database should have been seen completely
    std::sort(seenDatabases.begin(), seenDatabases.end());
    std::sort(expectedSeenDatabases.begin(), expectedSeenDatabases.end());
    EXPECT_EQ(seenDatabases, expectedSeenDatabases);
}
