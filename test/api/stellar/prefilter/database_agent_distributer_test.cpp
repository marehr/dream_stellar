#include <gtest/gtest.h>

#include <stellar/prefilter/database_agent_distributer.hpp>

TEST(DatabaseAgentDistributer, emptyDatabases)
{
    using TAlphabet = seqan::Dna5;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    std::vector<stellar::StellarDatabaseSegment<TAlphabet>> storage{};

    for (size_t splitCount = 1u; splitCount < 10u; ++splitCount)
    {
        auto segments = stellar::DatabaseAgentDistributer::distribute(storage, splitCount);
        EXPECT_EQ(segments.size(), 0u);
    }
}

TEST(DatabaseAgentDistributer, split9in3)
{
    using TAlphabet = seqan::Dna5;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGGGGGGGGG");

    std::vector<stellar::StellarDatabaseSegment<TAlphabet>> storage
    {
        {databases[0], 0u, 12u},
        {databases[0], 12u - 5u, 12u + 5u},
        {databases[0], 12u, 23u},
        {databases[1], 0u, 12u},
        {databases[1], 12u - 5u, 12u + 5u},
        {databases[1], 12u, 23u},
        {databases[2], 0u, 12u},
        {databases[2], 12u - 5u, 12u + 5u},
        {databases[2], 12u, 23u}
    };

    auto segments = stellar::DatabaseAgentDistributer::distribute(storage, 3u);

    // split into 3 groups a 3 segments
    ASSERT_EQ(segments.size(), 3u);
    EXPECT_EQ(segments[0].size(), 3u);
    EXPECT_EQ(segments[0][0], storage[0]);
    EXPECT_EQ(segments[0][1], storage[1]);
    EXPECT_EQ(segments[0][2], storage[2]);
    EXPECT_EQ(segments[1].size(), 3u);
    EXPECT_EQ(segments[1][0], storage[3]);
    EXPECT_EQ(segments[1][1], storage[4]);
    EXPECT_EQ(segments[1][2], storage[5]);
    EXPECT_EQ(segments[2].size(), 3u);
    EXPECT_EQ(segments[2][0], storage[6]);
    EXPECT_EQ(segments[2][1], storage[7]);
    EXPECT_EQ(segments[2][2], storage[8]);
}

TEST(DatabaseAgentDistributer, split9in4)
{
    using TAlphabet = seqan::Dna5;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGGGGGGGGG");

    std::vector<stellar::StellarDatabaseSegment<TAlphabet>> storage
    {
        {databases[0], 0u, 12u},
        {databases[0], 12u - 5u, 12u + 5u},
        {databases[0], 12u, 23u},
        {databases[1], 0u, 12u},
        {databases[1], 12u - 5u, 12u + 5u},
        {databases[1], 12u, 23u},
        {databases[2], 0u, 12u},
        {databases[2], 12u - 5u, 12u + 5u},
        {databases[2], 12u, 23u}
    };

    auto segments = stellar::DatabaseAgentDistributer::distribute(storage, 4u);

    // split into 4 groups a 2 segments (except first which has 3)
    ASSERT_EQ(segments.size(), 4u);
    EXPECT_EQ(segments[0].size(), 3u);
    EXPECT_EQ(segments[0][0], storage[0]);
    EXPECT_EQ(segments[0][1], storage[1]);
    EXPECT_EQ(segments[0][2], storage[2]);
    EXPECT_EQ(segments[1].size(), 2u);
    EXPECT_EQ(segments[1][0], storage[3]);
    EXPECT_EQ(segments[1][1], storage[4]);
    EXPECT_EQ(segments[2].size(), 2u);
    EXPECT_EQ(segments[2][0], storage[5]);
    EXPECT_EQ(segments[2][1], storage[6]);
    EXPECT_EQ(segments[3].size(), 2u);
    EXPECT_EQ(segments[3][0], storage[7]);
    EXPECT_EQ(segments[3][1], storage[8]);
}

TEST(DatabaseAgentDistributer, split9in5)
{
    using TAlphabet = seqan::Dna5;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGGGGGGGGG");

    std::vector<stellar::StellarDatabaseSegment<TAlphabet>> storage
    {
        {databases[0], 0u, 12u},
        {databases[0], 12u - 5u, 12u + 5u},
        {databases[0], 12u, 23u},
        {databases[1], 0u, 12u},
        {databases[1], 12u - 5u, 12u + 5u},
        {databases[1], 12u, 23u},
        {databases[2], 0u, 12u},
        {databases[2], 12u - 5u, 12u + 5u},
        {databases[2], 12u, 23u}
    };

    auto segments = stellar::DatabaseAgentDistributer::distribute(storage, 5u);

    // split into 5 groups a 1 segments (except first four which have 2 segments)
    ASSERT_EQ(segments.size(), 5u);
    EXPECT_EQ(segments[0].size(), 2u);
    EXPECT_EQ(segments[0][0], storage[0]);
    EXPECT_EQ(segments[0][1], storage[1]);
    EXPECT_EQ(segments[1].size(), 2u);
    EXPECT_EQ(segments[1][0], storage[2]);
    EXPECT_EQ(segments[1][1], storage[3]);
    EXPECT_EQ(segments[2].size(), 2u);
    EXPECT_EQ(segments[2][0], storage[4]);
    EXPECT_EQ(segments[2][1], storage[5]);
    EXPECT_EQ(segments[3].size(), 2u);
    EXPECT_EQ(segments[3][0], storage[6]);
    EXPECT_EQ(segments[3][1], storage[7]);
    EXPECT_EQ(segments[4].size(), 1u);
    EXPECT_EQ(segments[4][0], storage[8]);
}

TEST(DatabaseAgentDistributer, split9in9orMore)
{
    using TAlphabet = seqan::Dna5;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGGGGGGGGG");

    std::vector<stellar::StellarDatabaseSegment<TAlphabet>> storage
    {
        {databases[0], 0u, 12u},
        {databases[0], 12u - 5u, 12u + 5u},
        {databases[0], 12u, 23u},
        {databases[1], 0u, 12u},
        {databases[1], 12u - 5u, 12u + 5u},
        {databases[1], 12u, 23u},
        {databases[2], 0u, 12u},
        {databases[2], 12u - 5u, 12u + 5u},
        {databases[2], 12u, 23u}
    };

    for (size_t splitCount = 9u; splitCount < 20u; ++splitCount)
    {
        auto segments = stellar::DatabaseAgentDistributer::distribute(storage, splitCount);

        // split into 9 groups a 1 segments
        ASSERT_EQ(segments.size(), 9u);
        EXPECT_EQ(segments[0].size(), 1u);
        EXPECT_EQ(segments[0][0], storage[0]);
        EXPECT_EQ(segments[1].size(), 1u);
        EXPECT_EQ(segments[1][0], storage[1]);
        EXPECT_EQ(segments[2].size(), 1u);
        EXPECT_EQ(segments[2][0], storage[2]);
        EXPECT_EQ(segments[3].size(), 1u);
        EXPECT_EQ(segments[3][0], storage[3]);
        EXPECT_EQ(segments[4].size(), 1u);
        EXPECT_EQ(segments[4][0], storage[4]);
        EXPECT_EQ(segments[5].size(), 1u);
        EXPECT_EQ(segments[5][0], storage[5]);
        EXPECT_EQ(segments[6].size(), 1u);
        EXPECT_EQ(segments[6][0], storage[6]);
        EXPECT_EQ(segments[7].size(), 1u);
        EXPECT_EQ(segments[7][0], storage[7]);
        EXPECT_EQ(segments[8].size(), 1u);
        EXPECT_EQ(segments[8][0], storage[8]);
    }
}
