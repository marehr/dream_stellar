#include <gtest/gtest.h>

#include <stellar/test/alignment_fragment.hpp>

#include <stellar/test/fixture/512_simSeq_5e-2.hpp>

TEST(alignment, empty)
{
    using TAlphabet = seqan::Dna5;
    stellar::test::alignment_fragment<TAlphabet> fragment{"[]", "[]"};

    EXPECT_EQ(fragment.sequence1, "");
    EXPECT_EQ(fragment.sequence2, "");
    EXPECT_EQ(fragment.beginPosition1, 0);
    EXPECT_EQ(fragment.beginPosition2, 0);
    EXPECT_EQ(fragment.endPosition1, 0);
    EXPECT_EQ(fragment.endPosition2, 0);
    EXPECT_EQ(row(fragment.alignment, 0), "");
    EXPECT_EQ(row(fragment.alignment, 1), "");
}

TEST(alignment, empty2)
{
    using TAlphabet = seqan::Dna5;
    stellar::test::alignment_fragment<TAlphabet> fragment
    {
        "ACGTAC[]GT",
        "TGC[]ATGCA"
    };

    EXPECT_EQ(fragment.sequence1, "ACGTACGT");
    EXPECT_EQ(fragment.sequence2, "TGCATGCA");
    EXPECT_EQ(fragment.beginPosition1, 6);
    EXPECT_EQ(fragment.beginPosition2, 3);
    EXPECT_EQ(fragment.endPosition1, 6);
    EXPECT_EQ(fragment.endPosition2, 3);
    EXPECT_EQ(row(fragment.alignment, 0), "");
    EXPECT_EQ(row(fragment.alignment, 1), "");
}

TEST(alignment, singleChar)
{
    using TAlphabet = seqan::Dna5;
    stellar::test::alignment_fragment<TAlphabet> fragment
    {
        "ACGT[A]CGT",
        "TGC[A]TGCA"
    };

    EXPECT_EQ(fragment.sequence1, "ACGTACGT");
    EXPECT_EQ(fragment.sequence2, "TGCATGCA");
    EXPECT_EQ(fragment.beginPosition1, 4);
    EXPECT_EQ(fragment.beginPosition2, 3);
    EXPECT_EQ(fragment.endPosition1, 5);
    EXPECT_EQ(fragment.endPosition2, 4);
    EXPECT_EQ(row(fragment.alignment, 0), "A");
    EXPECT_EQ(row(fragment.alignment, 1), "A");
}

TEST(alignment, ungappedSequence)
{
    using TAlphabet = seqan::Dna5;
    stellar::test::alignment_fragment<TAlphabet> fragment
    {
        "ACGT[ACG]T",
        "TGC[ATG]CA"
    };

    EXPECT_EQ(fragment.sequence1, "ACGTACGT");
    EXPECT_EQ(fragment.sequence2, "TGCATGCA");
    EXPECT_EQ(fragment.beginPosition1, 4);
    EXPECT_EQ(fragment.beginPosition2, 3);
    EXPECT_EQ(fragment.endPosition1, 7);
    EXPECT_EQ(fragment.endPosition2, 6);
    EXPECT_EQ(row(fragment.alignment, 0), "ACG");
    EXPECT_EQ(row(fragment.alignment, 1), "ATG");
}

TEST(alignment, gappedSequence)
{
    using TAlphabet = seqan::Dna5;
    stellar::test::alignment_fragment<TAlphabet> fragment
    {
        "ACGT[-ACGT]",
        "TG""[CATG-]CA"
    };

    EXPECT_EQ(fragment.sequence1, "ACGTACGT");
    EXPECT_EQ(fragment.sequence2, "TGCATGCA");
    EXPECT_EQ(fragment.beginPosition1, 4);
    EXPECT_EQ(fragment.beginPosition2, 2);
    EXPECT_EQ(fragment.endPosition1, 8);
    EXPECT_EQ(fragment.endPosition2, 6);
    EXPECT_EQ(row(fragment.alignment, 0), "-ACGT");
    EXPECT_EQ(row(fragment.alignment, 1), "CATG-");
}

static stellar::test::alignment_fragment<seqan::Dna5> const globalFragment
{
      "[AAAAACCCCC]GGGGGTTTTT",
    "AA[--AAACCC--]CCGGGGGTTTTT",
};

TEST(alignment, gappedSequenceWithGlobalFragment)
{
    auto & fragment = globalFragment;

    EXPECT_EQ(fragment.sequence1, "AAAAACCCCCGGGGGTTTTT");
    EXPECT_EQ(fragment.sequence2, "AAAAACCCCCGGGGGTTTTT");
    EXPECT_EQ(fragment.beginPosition1, 0);
    EXPECT_EQ(fragment.beginPosition2, 2);
    EXPECT_EQ(fragment.endPosition1, 10);
    EXPECT_EQ(fragment.endPosition2, 8);
    EXPECT_EQ(row(fragment.alignment, 0), "AAAAACCCCC");
    EXPECT_EQ(row(fragment.alignment, 1), "--AAACCC--");
}

TEST(alignment, gappedSequenceWithWhitespace)
{
    using TAlphabet = seqan::Dna5;
    stellar::test::alignment_fragment<TAlphabet> fragment
    {
        "    [ AAAAA CCCCC ]  GGGGG  TTTTT  ",
        " A A[ --AAA CCC-- ]CCGGGGG  TTTTT ",
    };

    EXPECT_EQ(fragment.sequence1, "AAAAACCCCCGGGGGTTTTT");
    EXPECT_EQ(fragment.sequence2, "AAAAACCCCCGGGGGTTTTT");
    EXPECT_EQ(fragment.beginPosition1, 0);
    EXPECT_EQ(fragment.beginPosition2, 2);
    EXPECT_EQ(fragment.endPosition1, 10);
    EXPECT_EQ(fragment.endPosition2, 8);
    EXPECT_EQ(row(fragment.alignment, 0), "AAAAACCCCC");
    EXPECT_EQ(row(fragment.alignment, 1), "--AAACCC--");
}

TEST(alignment, gappedSequenceWithMultipleSequenceLines)
{
    using TAlphabet = seqan::Dna5;
    stellar::test::alignment_fragment<TAlphabet> const & fragment
        = stellar::test::fixture::a512_simSeq_0dot05::forward::fragment0;

    EXPECT_EQ(fragment.sequence1,
        "TTTGAAAAAACTTACTTGTCAATTT"
        "ATAAGAAGCGGGAATTTTGTATTGCTTCAATACAGTCTGTCGTTTGGGCACGACCAGTTCGGACCTCGCACCAATCGCTGACGATGCGGCAACGTACACA"
        "AACATGCGCATCTAACCGCTAGCTCACTATGACTATTTCGCTGATTATCTGGTGGGGGCATCCACTATCCGTATTGGCCGAAATACAGCGGCGTGATG"
        "AATCACGTCGGG"
        "AGTAGACGCCTCTCTTGTACTCGGG");

    EXPECT_EQ(fragment.sequence2,
        "AACGGGGTTGAACCAATTATTCATC"
        "ATAAGAAGCGGGAATTTTGTATTGCTACAATACAGTCTGTCTTTGGGCACGACCAGTTCGGACCTCGCACCAATCGCTGACGATGCGGCACGTACACA"
        "AACATGCGCATCTAACCGCTAGCTCACTATAGACTATTTCGCTGATTATTTGGTGGGGGCATCCACTATCCGTATTGGCCAAATACCAGCGGCGTGATG"
        "AATCCTCGAGG"
        "TTGTGGCATAAGCGAGCATACAGCG");
    EXPECT_EQ(fragment.beginPosition1, 25u);
    EXPECT_EQ(fragment.beginPosition2, 25u);
    EXPECT_EQ(fragment.endPosition1, 25u + /*alignment_length*/ 213 - /*gaps*/ 3u);
    EXPECT_EQ(fragment.endPosition2, 25u + /*alignment_length*/ 213 - /*gaps*/ 5u);
    EXPECT_EQ(row(fragment.alignment, 0),
        "ATAAGAAGCGGGAATTTTGTATTGCTTCAATACAGTCTGTCGTTTGGGCACGACCAGTTCGGACCTCGCACCAATCGCTGACGATGCGGCAACGTACACA"
        "AACATGCGCATCTAACCGCTAGCTCACTAT-GACTATTTCGCTGATTATCTGGTGGGGGCATCCACTATCCGTATTGGCCGAAATA-CAGCGGCGTGATG"
        "AATCACGTCG-GG");
    EXPECT_EQ(row(fragment.alignment, 1),
        "ATAAGAAGCGGGAATTTTGTATTGCTACAATACAGTCTGTC-TTTGGGCACGACCAGTTCGGACCTCGCACCAATCGCTGACGATGCGGC-ACGTACACA"
        "AACATGCGCATCTAACCGCTAGCTCACTATAGACTATTTCGCTGATTATTTGGTGGGGGCATCCACTATCCGTATTGGCC-AAATACCAGCGGCGTGATG"
        "AATC-C-TCGAGG");
}
