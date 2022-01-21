#include <gtest/gtest.h>

#include <stellar/prefilter/nsegment_database_agent_splitter.hpp>

#include <stellar/test/prefilter/expect_storage.hpp>

TEST(NSegmentDatabaseAgentSplitter, constructor)
{
    {
        stellar::NSegmentDatabaseAgentSplitter const splitter{};
        EXPECT_EQ(splitter.segmentCount, 1u);
    }

    {
        stellar::NSegmentDatabaseAgentSplitter const splitter{18u};
        EXPECT_EQ(splitter.segmentCount, 18u);
    }
}

TEST(NSegmentDatabaseAgentSplitter, emptyDatabases)
{
    using TAlphabet = seqan::Dna5;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};

    for (size_t minLength = 0u; minLength < 10u; ++minLength)
    {
        stellar::NSegmentDatabaseAgentSplitter const splitter{2u};
        auto storage = splitter.split(databases, minLength);
        EXPECT_EQ(storage.size(), 0u);
        EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));
    }
}

TEST(NSegmentDatabaseAgentSplitter, divideDatabaseIn1SegmentWithMinLengthX)
{
    using TAlphabet = seqan::Dna5;
    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGG");

    EXPECT_EQ(length(databases[0]), 23u);
    EXPECT_EQ(length(databases[1]), 19u);
    EXPECT_EQ(length(databases[2]), 16u);

    stellar::NSegmentDatabaseAgentSplitter const splitter{1u};

    for (size_t minLength = 0u; minLength < 24u; ++minLength)
    {
        auto storage = splitter.split(databases, minLength);

        // all most 3 segmented sequences
        ASSERT_EQ(storage.size(), (minLength <= 16u) + (minLength <= 19u) + (minLength <= 23u));
        EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));

        if (storage.size() >= 1u) {
            EXPECT_EQ(storage[0], (TDatabaseSegment{databases[0], 0u, 23u}));
        }
        if (storage.size() >= 2u) {
            EXPECT_EQ(storage[1], (TDatabaseSegment{databases[1], 0u, 19u}));
        }
        if (storage.size() >= 3u) {
            EXPECT_EQ(storage[2], (TDatabaseSegment{databases[2], 0u, 16u}));
        }
    }
}

TEST(NSegmentDatabaseAgentSplitter, divideDatabaseIn2SegmentsWithMinLength5)
{
    using TAlphabet = seqan::Dna5;
    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGG");

    EXPECT_EQ(length(databases[0]), 23u);
    EXPECT_EQ(length(databases[1]), 19u);
    EXPECT_EQ(length(databases[2]), 16u);

    size_t const minLength = 5u;
    stellar::NSegmentDatabaseAgentSplitter const splitter{2u};
    auto storage = splitter.split(databases, minLength);

    // all together 9 segmented sequences
    ASSERT_EQ(storage.size(), 9u);
    EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));

    size_t const segmentLength0 = 23u / 2u; // remainder 1u
    EXPECT_EQ(segmentLength0, 11u);
    EXPECT_EQ(storage[0], (TDatabaseSegment{databases[0], 0u, 12u})); // length: 12
    EXPECT_EQ(storage[1], (TDatabaseSegment{databases[0], 12u - minLength, 12u + minLength}));
    EXPECT_EQ(storage[2], (TDatabaseSegment{databases[0], 12u, 23u})); // length: 11

    size_t const segmentLength1 = 19u / 2u; // remainder 1u
    EXPECT_EQ(segmentLength1, 9u);
    EXPECT_EQ(storage[3], (TDatabaseSegment{databases[1], 0u, 10u})); // length: 10
    EXPECT_EQ(storage[4], (TDatabaseSegment{databases[1], 10u - minLength, 10u + minLength}));
    EXPECT_EQ(storage[5], (TDatabaseSegment{databases[1], 10u, 19u})); // length: 9

    size_t const segmentLength2 = 16u / 2u; // no remainder
    EXPECT_EQ(segmentLength2, 8u);
    EXPECT_EQ(storage[6], (TDatabaseSegment{databases[2], 0u, 8u})); // length: 8
    EXPECT_EQ(storage[7], (TDatabaseSegment{databases[2], 8u - minLength, 8u + minLength}));
    EXPECT_EQ(storage[8], (TDatabaseSegment{databases[2], 8u, 16u})); // length: 8
}

TEST(NSegmentDatabaseAgentSplitter, divideDatabaseIn2SegmentsWithMinLength0)
{
    using TAlphabet = seqan::Dna5;
    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGG");

    EXPECT_EQ(length(databases[0]), 23u);
    EXPECT_EQ(length(databases[1]), 19u);
    EXPECT_EQ(length(databases[2]), 16u);

    size_t const minLength = 0u;
    stellar::NSegmentDatabaseAgentSplitter const splitter{2u};
    auto storage = splitter.split(databases, minLength);

    // all together 6 segmented sequences
    ASSERT_EQ(storage.size(), 6u);
    EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));

    size_t const segmentLength0 = 23u / 2u; // remainder 1u
    EXPECT_EQ(segmentLength0, 11u);
    EXPECT_EQ(storage[0], (TDatabaseSegment{databases[0], 0u, 12u})); // length: 12
    EXPECT_EQ(storage[1], (TDatabaseSegment{databases[0], 12u, 23u})); // length: 11

    size_t const segmentLength1 = 19u / 2u; // remainder 1u
    EXPECT_EQ(segmentLength1, 9u);
    EXPECT_EQ(storage[2], (TDatabaseSegment{databases[1], 0u, 10u})); // length: 10
    EXPECT_EQ(storage[3], (TDatabaseSegment{databases[1], 10u, 19u})); // length: 9

    size_t const segmentLength2 = 16u / 2u; // no remainder
    EXPECT_EQ(segmentLength2, 8u);
    EXPECT_EQ(storage[4], (TDatabaseSegment{databases[2], 0u, 8u})); // length: 8
    EXPECT_EQ(storage[5], (TDatabaseSegment{databases[2], 8u, 16u})); // length: 8
}

TEST(NSegmentDatabaseAgentSplitter, divideDatabaseIn3SegmentsWithMinLength5)
{
    using TAlphabet = seqan::Dna5;
    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGG");

    EXPECT_EQ(length(databases[0]), 23u);
    EXPECT_EQ(length(databases[1]), 19u);
    EXPECT_EQ(length(databases[2]), 16u);

    size_t const minLength = 5u;
    stellar::NSegmentDatabaseAgentSplitter const splitter{3u};
    auto storage = splitter.split(databases, minLength);

    // all together 15 segmented sequences
    ASSERT_EQ(storage.size(), 15u);
    EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));

    size_t const segmentLength0 = 23u / 3u; // remainder 2
    EXPECT_EQ(segmentLength0, 7u);
    EXPECT_EQ(storage[0], (TDatabaseSegment{databases[0], 0u, 8u})); // length: 8
    EXPECT_EQ(storage[1], (TDatabaseSegment{databases[0], 8u - minLength, 8u + minLength}));
    EXPECT_EQ(storage[2], (TDatabaseSegment{databases[0], 8u, 16u})); // length: 8
    EXPECT_EQ(storage[3], (TDatabaseSegment{databases[0], 16u - minLength, 16u + minLength}));
    EXPECT_EQ(storage[4], (TDatabaseSegment{databases[0], 16u, 23u})); // length: 7

    size_t const segmentLength1 = 19u / 3u; // remainder 1
    EXPECT_EQ(segmentLength1, 6u);
    EXPECT_EQ(storage[5], (TDatabaseSegment{databases[1], 0u, 7u})); // length: 7
    EXPECT_EQ(storage[6], (TDatabaseSegment{databases[1], 7u - 5u, 7u + 5u}));
    EXPECT_EQ(storage[7], (TDatabaseSegment{databases[1], 7u, 13u})); // length: 6
    EXPECT_EQ(storage[8], (TDatabaseSegment{databases[1], 13u - 5u, 13u + 5u}));
    EXPECT_EQ(storage[9], (TDatabaseSegment{databases[1], 13u, 19u})); // length: 6

    size_t const segmentLength2 = 16u / 3u; // remainder 1
    EXPECT_EQ(segmentLength2, 5u);
    EXPECT_EQ(storage[10], (TDatabaseSegment{databases[2], 0u, 6u})); // length: 6
    EXPECT_EQ(storage[11], (TDatabaseSegment{databases[2], 6u - 5u, 6u + 5u}));
    EXPECT_EQ(storage[12], (TDatabaseSegment{databases[2], 6u, 11u})); // length: 5
    EXPECT_EQ(storage[13], (TDatabaseSegment{databases[2], 11u - 5u, 11u + 5u}));
    EXPECT_EQ(storage[14], (TDatabaseSegment{databases[2], 11u, 16u})); // length: 5
}

TEST(NSegmentDatabaseAgentSplitter, divideDatabaseIn4SegmentsWithMinLength5)
{
    using TAlphabet = seqan::Dna5;
    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGG");

    EXPECT_EQ(length(databases[0]), 23u);
    EXPECT_EQ(length(databases[1]), 19u);
    EXPECT_EQ(length(databases[2]), 16u);

    size_t const minLength = 5u;
    stellar::NSegmentDatabaseAgentSplitter const splitter{4u};
    auto storage = splitter.split(databases, minLength);

    // all together 17 segmented sequences
    ASSERT_EQ(storage.size(), 17u);
    EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));

    size_t const segmentLength0 = 23u / 4u; // remainder 3
    EXPECT_EQ(segmentLength0, 5u);
    EXPECT_EQ(storage[0], (TDatabaseSegment{databases[0], 0u, 6u})); // length: 6
    EXPECT_EQ(storage[1], (TDatabaseSegment{databases[0], 6u - minLength, 6u + minLength}));
    EXPECT_EQ(storage[2], (TDatabaseSegment{databases[0], 6u, 12u})); // length: 6
    EXPECT_EQ(storage[3], (TDatabaseSegment{databases[0], 12u - minLength, 12u + minLength}));
    EXPECT_EQ(storage[4], (TDatabaseSegment{databases[0], 12u, 18u})); // length: 6
    EXPECT_EQ(storage[5], (TDatabaseSegment{databases[0], 18u - minLength, 18u + minLength}));
    EXPECT_EQ(storage[6], (TDatabaseSegment{databases[0], 18u, 23u})); // length: 5

    // edge case: segmentLength1 is below minLength, set segmentLength to at least minLength
    size_t const segmentLength1 = 19u / 4u; // remainder 3
    EXPECT_EQ(segmentLength1, 4u);
    size_t const newSegmentCount1 = 19u / minLength;
    EXPECT_EQ(newSegmentCount1, 3);
    size_t const newSegmentLength1 = 19u / newSegmentCount1; // remainder 1
    EXPECT_EQ(newSegmentLength1, 6u);
    EXPECT_EQ(storage[7+0], (TDatabaseSegment{databases[1], 0u, 7u})); // length: 7
    EXPECT_EQ(storage[7+1], (TDatabaseSegment{databases[1], 7u - minLength, 7u + minLength}));
    EXPECT_EQ(storage[7+2], (TDatabaseSegment{databases[1], 7u, 13u})); // length: 6
    EXPECT_EQ(storage[7+3], (TDatabaseSegment{databases[1], 13u - minLength, 13u + minLength}));
    EXPECT_EQ(storage[7+4], (TDatabaseSegment{databases[1], 13u, 19u})); // length: 6

    // edge case: segmentLength2 is below minLength, set segmentLength to at least minLength
    size_t const segmentLength2 = 16u / 4u; // no remainder
    EXPECT_EQ(segmentLength2, 4u);
    size_t const newSegmentCount2 = 16u / minLength;
    EXPECT_EQ(newSegmentCount2, 3);
    size_t const newSegmentLength2 = 16u / newSegmentCount2; // remainder 1
    EXPECT_EQ(newSegmentLength2, 5u);
    EXPECT_EQ(storage[7+5+0], (TDatabaseSegment{databases[2], 0u, 6u})); // length: 6
    EXPECT_EQ(storage[7+5+1], (TDatabaseSegment{databases[2], 6u - minLength, 6u + minLength}));
    EXPECT_EQ(storage[7+5+2], (TDatabaseSegment{databases[2], 6u, 11u})); // length: 5
    EXPECT_EQ(storage[7+5+3], (TDatabaseSegment{databases[2], 11u - minLength, 11u + minLength}));
    EXPECT_EQ(storage[7+5+4], (TDatabaseSegment{databases[2], 11u, 16u})); // length: 5
}

TEST(NSegmentDatabaseAgentSplitter, divideDatabaseIn4SegmentsWithMinLength10)
{
    using TAlphabet = seqan::Dna5;
    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGG");

    EXPECT_EQ(length(databases[0]), 23u);
    EXPECT_EQ(length(databases[1]), 19u);
    EXPECT_EQ(length(databases[2]), 16u);

    size_t const minLength = 10u;
    stellar::NSegmentDatabaseAgentSplitter const splitter{4u};
    auto storage = splitter.split(databases, minLength);

    // all together 5 segmented sequences
    ASSERT_EQ(storage.size(), 5u);
    EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));

    // edge case: segmentLength0 is below minLength, set segmentLength to at least minLength
    size_t const segmentLength0 = 23u / 4u; // remainder 3
    EXPECT_EQ(segmentLength0, 5u);
    size_t const newSegmentCount0 = 23u / minLength;
    EXPECT_EQ(newSegmentCount0, 2);
    size_t const newSegmentLength0 = 23u / newSegmentCount0; // remainder 1
    EXPECT_EQ(newSegmentLength0, 11u);
    EXPECT_EQ(storage[0], (TDatabaseSegment{databases[0], 0u, 12u})); // length: 12
    EXPECT_EQ(storage[1], (TDatabaseSegment{databases[0], 12u - minLength, 12u + minLength}));
    EXPECT_EQ(storage[2], (TDatabaseSegment{databases[0], 12u, 23u})); // length: 11

    // edge case: segmentLength1 is below minLength, set segmentLength to at least minLength
    size_t const segmentLength1 = 19u / 4u; // remainder 3
    EXPECT_EQ(segmentLength1, 4u);
    size_t const newSegmentCount1 = 19u / minLength;
    EXPECT_EQ(newSegmentCount1, 1);
    size_t const newSegmentLength1 = 19u / newSegmentCount1; // remainder 0
    EXPECT_EQ(newSegmentLength1, 19u);
    EXPECT_EQ(storage[3], (TDatabaseSegment{databases[1], 0u, 19u})); // length: 19

    // edge case: segmentLength1 is below minLength, set segmentLength to at least minLength
    size_t const segmentLength2 = 16u / 4u; // no remainder
    EXPECT_EQ(segmentLength2, 4u);
    size_t const newSegmentCount2 = 16u / minLength;
    EXPECT_EQ(newSegmentCount2, 1);
    size_t const newSegmentLength2 = 16u / newSegmentCount2; // remainder 0
    EXPECT_EQ(newSegmentLength2, 16u);
    EXPECT_EQ(storage[4], (TDatabaseSegment{databases[2], 0u, 16u})); // length: 16
}

TEST(NSegmentDatabaseAgentSplitter, divideDatabaseIn4SegmentsWithMinLength19)
{
    using TAlphabet = seqan::Dna5;
    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGG");

    EXPECT_EQ(length(databases[0]), 23u);
    EXPECT_EQ(length(databases[1]), 19u);
    EXPECT_EQ(length(databases[2]), 16u);

    size_t const minLength = 19u;
    stellar::NSegmentDatabaseAgentSplitter const splitter{4u};
    auto storage = splitter.split(databases, minLength);

    // all together 2 segmented sequences
    ASSERT_EQ(storage.size(), 2u);
    EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));

    // edge case: segmentLength0 is below minLength, set segmentLength to at least minLength
    size_t const segmentLength0 = 23u / 4u; // remainder 3
    EXPECT_EQ(segmentLength0, 5u);
    size_t const newSegmentCount0 = 23u / minLength;
    EXPECT_EQ(newSegmentCount0, 1);
    size_t const newSegmentLength0 = 23u / newSegmentCount0; // remainder 0
    EXPECT_EQ(newSegmentLength0, 23u);
    EXPECT_EQ(storage[0], (TDatabaseSegment{databases[0], 0u, 23u})); // length: 23

    // edge case: database length is exactly minLength
    size_t const segmentLength1 = 19u / 4u; // remainder 3
    EXPECT_EQ(segmentLength1, 4u);
    size_t const newSegmentCount1 = 19u / minLength;
    EXPECT_EQ(newSegmentCount1, 1);
    size_t const newSegmentLength1 = 19u / newSegmentCount1; // remainder 0
    EXPECT_EQ(newSegmentLength1, 19u);
    EXPECT_EQ(storage[1], (TDatabaseSegment{databases[1], 0u, 19u})); // length: 19

    // edge case: database is shorter than minLength
    size_t const segmentLength2 = 16u / 4u; // no remainder
    EXPECT_EQ(segmentLength2, 4u);
    size_t const newSegmentCount2 = 16u / minLength;
    EXPECT_EQ(newSegmentCount2, 0); // this sequence does not fulfill min length requirement
}
