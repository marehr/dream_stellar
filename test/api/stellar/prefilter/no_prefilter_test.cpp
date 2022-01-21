#include <gtest/gtest.h>

#include <stellar/prefilter/no_query_prefilter.hpp>
#include <stellar/prefilter/whole_database_agent_splitter.hpp>
#include <stellar/prefilter/nsegment_database_agent_splitter.hpp>

TEST(NoQueryPrefilter, WholeDatabaseAgentSplitter)
{
    using TAlphabet = seqan::Dna5;
    using TQueryFilter = seqan::StringSet<seqan::String<TAlphabet> > const &;
    using TAgentSplitter = stellar::WholeDatabaseAgentSplitter;

    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGGGGGGGGG");

    seqan::StringSet<seqan::String<TAlphabet> > queries{};
    appendValue(queries, "AAAA");
    appendValue(queries, "CCCC");
    appendValue(queries, "GGGG");
    appendValue(queries, "TTTT");

    // Note: TQueryFilter are the queries themselve to check whith which queries the actual filter would have been
    //       constructed.
    using TPrefilter = stellar::NoQueryPrefilter<TAlphabet, TQueryFilter, TAgentSplitter>;
    using TPrefilterAgent = TPrefilter::Agent;

    size_t const agentCount = 2u;
    size_t const minLength = 2u;
    TAgentSplitter const splitter{};
    TPrefilter prefilter{databases, queries};

    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    std::vector<TDatabaseSegment> seenDatabases{};
    std::vector<TDatabaseSegment> expectedSeenDatabases = splitter.split(databases, minLength);

    for (TPrefilterAgent prefilterAgent: prefilter.agents(agentCount, minLength))
    {
        using TDatabaseSegments = typename TPrefilter::TDatabaseSegments;
        using TQueryFilter = typename TPrefilter::TQueryFilter;

        prefilterAgent.prefilter([&](TDatabaseSegments databaseSegments, TQueryFilter queryFilter)
        {
            // collect all seen databases
            for (TDatabaseSegment const & databaseSegment : databaseSegments)
                seenDatabases.push_back(databaseSegment);

            // This ensures that the queryFilter would have been constructed with ALL queries
            EXPECT_EQ(queryFilter, queries);
        });
    }

    // Note: every database should have been seen completely
    std::sort(seenDatabases.begin(), seenDatabases.end());
    std::sort(expectedSeenDatabases.begin(), expectedSeenDatabases.end());
    EXPECT_EQ(seenDatabases, expectedSeenDatabases);
}

TEST(NoQueryPrefilter, NSegmentDatabaseAgentSplitter)
{
    using TAlphabet = seqan::Dna5;
    using TQueryFilter = seqan::StringSet<seqan::String<TAlphabet> > const &;
    using TAgentSplitter = stellar::NSegmentDatabaseAgentSplitter;

    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGGGGGGGGG");

    seqan::StringSet<seqan::String<TAlphabet> > queries{};
    appendValue(queries, "AAAA");
    appendValue(queries, "CCCC");
    appendValue(queries, "GGGG");
    appendValue(queries, "TTTT");

    // Note: TQueryFilter are the queries themselve to check whith which queries the actual filter would have been
    //       constructed.
    using TPrefilter = stellar::NoQueryPrefilter<TAlphabet, TQueryFilter, TAgentSplitter>;
    using TPrefilterAgent = TPrefilter::Agent;

    size_t const agentCount = 3u;
    size_t const minLength = 2u;
    TAgentSplitter const splitter{4u};
    TPrefilter prefilter{databases, queries, splitter};

    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    std::vector<TDatabaseSegment> seenDatabases{};
    std::vector<TDatabaseSegment> expectedSeenDatabases = splitter.split(databases, minLength);

    for (TPrefilterAgent prefilterAgent: prefilter.agents(agentCount, minLength))
    {
        using TDatabaseSegments = typename TPrefilter::TDatabaseSegments;
        using TQueryFilter = typename TPrefilter::TQueryFilter;

        prefilterAgent.prefilter([&](TDatabaseSegments databaseSegments, TQueryFilter queryFilter)
        {
            // collect all seen databases
            for (TDatabaseSegment const & databaseSegment : databaseSegments)
                seenDatabases.push_back(databaseSegment);

            // This ensures that the queryFilter would have been constructed with ALL queries
            EXPECT_EQ(queryFilter, queries);
        });
    }

    // Note: every database should have been seen completely
    std::sort(seenDatabases.begin(), seenDatabases.end());
    std::sort(expectedSeenDatabases.begin(), expectedSeenDatabases.end());
    EXPECT_EQ(seenDatabases, expectedSeenDatabases);
}
