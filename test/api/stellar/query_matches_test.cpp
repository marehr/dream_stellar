#include <gtest/gtest.h>

#include <string_view>

#include <stellar/stellar.hpp>

#include <stellar/test/alignment_fragment.hpp>

using TAlphabet = seqan::Dna5;
using TStellarMatch = stellar::StellarMatch<seqan::String<TAlphabet> const, seqan::CharString>;

TEST(StellarMatch, types)
{
    static_assert(std::is_same<typename TStellarMatch::TRow,
                               seqan::Gaps<seqan::String<TAlphabet> const, seqan::ArrayGaps>>::value, "");
}

TEST(maskOverlaps, empty)
{
    stellar::QueryMatches<TStellarMatch> matches{};

    EXPECT_EQ(length(matches.matches), 0u);
    EXPECT_EQ(matches.disabled, false);
    EXPECT_EQ(matches.lengthAdjustment, 0u);

    maskOverlaps(matches.matches, 10u);

    EXPECT_EQ(length(matches.matches), 0u);
    EXPECT_EQ(matches.disabled, false);
    EXPECT_EQ(matches.lengthAdjustment, 0u);
}

// Note: maskOverlaps assumes that all matches have at least the same
// TEST(maskOverlaps, minLengthNotSatisfied)
// {
//     TSequence database{"AAAAACCCCCGGGGGTTTTT"};
//     TSequence query{"AAAAA"};
//
//     stellar::QueryMatches<TStellarMatch> matches{};
//
//     {
//         TAlign alignment = make_alignment(seqan::infix(database, 0, 5),
//                                           seqan::infix(query, 0, 5),
//                                           "AAAAA",
//                                           "AAAAA");
//
//         TStellarMatch const match{alignment, "database1", false};
//         appendValue(matches.matches, match);
//     }
//
//     // min length 15 won't be satisfied by any match
//     maskOverlaps(matches.matches, 15u);
//
//     ASSERT_EQ(length(matches.matches), 1u);
//     EXPECT_EQ(matches.matches[0].id, "database1");
//     EXPECT_EQ(matches.matches[0].orientation, false);
//     EXPECT_EQ(matches.matches[0].begin1, 0);
//     EXPECT_EQ(matches.matches[0].end1, 5);
//     EXPECT_EQ(matches.matches[0].row1, "AAAAA");
//     EXPECT_EQ(matches.matches[0].begin2, 0);
//     EXPECT_EQ(matches.matches[0].end2, 5);
//     EXPECT_EQ(matches.matches[0].row2, "AAAAA");
//
//     EXPECT_EQ(matches.disabled, false);
// }

template <typename TAlphabet, typename TId>
void expect_eq_fixture(stellar::StellarMatch<seqan::String<TAlphabet> const, TId> const & match,
                       stellar::test::alignment_fragment<TAlphabet> const & fixture)
{
    EXPECT_EQ(match.begin1, fixture.beginPosition1);
    EXPECT_EQ(match.end1, fixture.endPosition1);
    EXPECT_EQ(match.row1, row(fixture.alignment, 0u));
    EXPECT_EQ(match.begin2, fixture.beginPosition2);
    EXPECT_EQ(match.end2, fixture.endPosition2);
    EXPECT_EQ(match.row2, row(fixture.alignment, 1u));
}

TEST(maskOverlaps, duplicates)
{
    std::vector<stellar::test::alignment_fragment<TAlphabet>> fixtures;
    stellar::test::alignment_fragment<TAlphabet> fixture1
    {
        "[AAAAA]CCCCCGGGGGTTTTT",
        "[AAAAA]",
    };
    fixtures.push_back(fixture1);
    fixtures.push_back(fixture1);

    stellar::QueryMatches<TStellarMatch> matches{};

    for (stellar::test::alignment_fragment<TAlphabet> const & fixture : fixtures)
    {
        TStellarMatch const match{fixture.alignment, "database1", false};
        appendValue(matches.matches, match);
    }

    maskOverlaps(matches.matches, 5u);

    ASSERT_EQ(length(matches.matches), 2u);
    EXPECT_EQ(matches.matches[0].id, "database1");
    EXPECT_EQ(matches.matches[0].orientation, false);
    expect_eq_fixture(matches.matches[0], fixtures[0]);

    EXPECT_EQ(matches.matches[1].id, TStellarMatch::INVALID_ID);
    EXPECT_EQ(matches.matches[1].orientation, false);
    expect_eq_fixture(matches.matches[1], fixtures[1]);

    EXPECT_EQ(matches.disabled, false);
}

TEST(maskOverlaps, nonOverlap)
{
    std::vector<stellar::test::alignment_fragment<TAlphabet>> fixtures;
    stellar::test::alignment_fragment<TAlphabet> fixture1
    {
        "[AAAAA]CCCCCGGGGGTTTTT",
        "[AAAAA]",
    };
    stellar::test::alignment_fragment<TAlphabet> fixture2
    {
        "AAAAA[CCCCC]GGGGGTTTTT",
        "[AAAAA]"
    };
    fixtures.push_back(fixture1);
    fixtures.push_back(fixture2);

    stellar::QueryMatches<TStellarMatch> matches{};

    for (stellar::test::alignment_fragment<TAlphabet> const & fixture : fixtures)
    {
        TStellarMatch const match{fixture.alignment, "database1", false};
        appendValue(matches.matches, match);
    }

    maskOverlaps(matches.matches, 5u);

    ASSERT_EQ(length(matches.matches), 2u);
    EXPECT_EQ(matches.matches[0].id, "database1");
    EXPECT_EQ(matches.matches[0].orientation, false);
    expect_eq_fixture(matches.matches[0], fixtures[0]);

    EXPECT_EQ(matches.matches[1].id, "database1");
    EXPECT_EQ(matches.matches[1].orientation, false);
    expect_eq_fixture(matches.matches[1], fixtures[1]);

    EXPECT_EQ(matches.disabled, false);
}

stellar::test::alignment_fragment<TAlphabet> containmentFixture1
{
      "[AAAAACCCCC]GGGGGTTTTT",
    "AA[--AAACCC--]CCGGGGGTTTTT",
};

stellar::test::alignment_fragment<TAlphabet> containmentFixture2
{
      "A[AAAACCCC]CGGGGGTTTTT",
    "AAA[--AACC--]CCCGGGGGTTTTT",
};

TEST(maskOverlaps, gapped_sequence3_contained_in_gapped_sequence1_and_gapped_sequence4_contained_in_gapped_sequence2)
{
    stellar::QueryMatches<TStellarMatch> matches{};
    std::vector<stellar::test::alignment_fragment<TAlphabet>> fixtures{containmentFixture1, containmentFixture2};

    for (stellar::test::alignment_fragment<TAlphabet> const & fixture : fixtures)
    {
        TStellarMatch const match{fixture.alignment, "database1", false};
        appendValue(matches.matches, match);
    }

    maskOverlaps(matches.matches, 5u);

    ASSERT_EQ(length(matches.matches), 2u);
    EXPECT_EQ(matches.matches[0].id, "database1");
    EXPECT_EQ(matches.matches[0].orientation, false);
    expect_eq_fixture(matches.matches[0], fixtures[0]);

    EXPECT_EQ(matches.matches[1].id, TStellarMatch::INVALID_ID);
    EXPECT_EQ(matches.matches[1].orientation, false);
    expect_eq_fixture(matches.matches[1], fixtures[1]);

    EXPECT_EQ(matches.disabled, false);
}

TEST(maskOverlaps, gapped_sequence1_contained_in_gapped_sequence3_and_gapped_sequence2_contained_in_gapped_sequence4)
{
    stellar::QueryMatches<TStellarMatch> matches{};
    std::vector<stellar::test::alignment_fragment<TAlphabet>> fixtures{containmentFixture2, containmentFixture1};

    for (stellar::test::alignment_fragment<TAlphabet> const & fixture : fixtures)
    {
        TStellarMatch const match{fixture.alignment, "database1", false};
        appendValue(matches.matches, match);
    }

    maskOverlaps(matches.matches, 5u);

    ASSERT_EQ(length(matches.matches), 2u);
    // note: matches will be sorted
    EXPECT_EQ(matches.matches[0].id, "database1");
    EXPECT_EQ(matches.matches[0].orientation, false);
    expect_eq_fixture(matches.matches[0], fixtures[1]);

    EXPECT_EQ(matches.matches[1].id, TStellarMatch::INVALID_ID);
    EXPECT_EQ(matches.matches[1].orientation, false);
    expect_eq_fixture(matches.matches[1], fixtures[0]);

    EXPECT_EQ(matches.disabled, false);
}
