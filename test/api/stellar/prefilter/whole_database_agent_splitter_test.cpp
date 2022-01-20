#include <gtest/gtest.h>

#include <stellar/prefilter/whole_database_agent_splitter.hpp>

#include <stellar/test/prefilter/expect_storage.hpp>

TEST(WholeDatabaseAgentSplitter, emptyDatabases)
{
    using TAlphabet = seqan::Dna5;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};

    for (size_t minLength = 0u; minLength < 10u; ++minLength)
    {
        stellar::WholeDatabaseAgentSplitter const splitter{};
        auto storage = splitter.split(databases, minLength);
        EXPECT_EQ(storage.size(), 0u);
        EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));
    }
}

TEST(WholeDatabaseAgentSplitter, databases2)
{
    using TAlphabet = seqan::Dna5;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCCCCCC");

    size_t const minLength = 5u;
    stellar::WholeDatabaseAgentSplitter const splitter{};
    auto storage = splitter.split(databases, minLength);
    // all together 2 sequences
    ASSERT_EQ(storage.size(), 2u);
    EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));

    EXPECT_EQ(storage[0].asInfixSegment(), databases[0]);
    EXPECT_EQ(storage[1].asInfixSegment(), databases[1]);
}

TEST(WholeDatabaseAgentSplitter, databases3)
{
    using TAlphabet = seqan::Dna5;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGGGGGGGGG");

    size_t const minLength = 5u;
    stellar::WholeDatabaseAgentSplitter const splitter{};
    auto storage = splitter.split(databases, minLength);
    // all together 3 sequences
    ASSERT_EQ(storage.size(), 3u);
    EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));

    EXPECT_EQ(storage[0].asInfixSegment(), databases[0]);
    EXPECT_EQ(storage[1].asInfixSegment(), databases[1]);
    EXPECT_EQ(storage[2].asInfixSegment(), databases[2]);
}

TEST(WholeDatabaseAgentSplitter, databases4)
{
    using TAlphabet = seqan::Dna5;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGGGGGGGGGGGGGGG");
    appendValue(databases, "TTTTTTTTTTTTTTTTTTTTTTT");

    size_t const minLength = 5u;
    stellar::WholeDatabaseAgentSplitter const splitter{};
    auto storage = splitter.split(databases, minLength);
    // all together 4 sequences
    ASSERT_EQ(storage.size(), 4u);
    EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));

    EXPECT_EQ(storage[0].asInfixSegment(), databases[0]);
    EXPECT_EQ(storage[1].asInfixSegment(), databases[1]);
    EXPECT_EQ(storage[2].asInfixSegment(), databases[2]);
    EXPECT_EQ(storage[3].asInfixSegment(), databases[3]);
}

TEST(WholeDatabaseAgentSplitter, databases5)
{
    using TAlphabet = seqan::Dna5;
    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "AAAAAAAAAAAAAAAAAAAAAAA");
    appendValue(databases, "CCCCCCCCCCCCCCC");
    appendValue(databases, "GGGGGGGGGG");
    appendValue(databases, "TTTTTTTTTTTTTTTTT");
    appendValue(databases, "ACGTACGTA");

    EXPECT_EQ(length(databases[0]), 23u); // >= 15
    EXPECT_EQ(length(databases[1]), 15u); // >= 15
    EXPECT_EQ(length(databases[2]), 10u); // < 15
    EXPECT_EQ(length(databases[3]), 17u); // >= 15
    EXPECT_EQ(length(databases[4]), 9u); // < 15

    size_t const minLength = 15u;
    stellar::WholeDatabaseAgentSplitter const splitter{};
    auto storage = splitter.split(databases, minLength);
    // all together 3 sequences
    ASSERT_EQ(storage.size(), 3u);
    EXPECT_TRUE(stellar::test::expectStorage(storage, minLength));

    EXPECT_EQ(storage[0].asInfixSegment(), databases[0]);
    EXPECT_EQ(storage[1].asInfixSegment(), databases[1]);
    EXPECT_EQ(storage[2].asInfixSegment(), databases[3]);
}
