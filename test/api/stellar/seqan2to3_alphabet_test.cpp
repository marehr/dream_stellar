#include <gtest/gtest.h>

#include <stellar/seqan2to3_alphabet.hpp>

TEST(dna5, concepts)
{
    EXPECT_TRUE((seqan3::semialphabet<seqan::Dna5>));
    EXPECT_TRUE((seqan3::alphabet<seqan::Dna5>));
    EXPECT_TRUE((seqan3::nucleotide_alphabet<seqan::Dna5>));
}

TEST(dna5, to_rank)
{
    seqan::Dna5 dna5[] = {'A', 'C', 'G', 'T', 'N'};

    EXPECT_EQ(seqan3::to_rank(dna5[0]), 0);
    EXPECT_EQ(seqan3::to_rank(dna5[1]), 1);
    EXPECT_EQ(seqan3::to_rank(dna5[2]), 2);
    EXPECT_EQ(seqan3::to_rank(dna5[3]), 3);
    EXPECT_EQ(seqan3::to_rank(dna5[4]), 4);
}

TEST(dna5, to_char)
{
    seqan::Dna5 dna5[] = {'A', 'C', 'G', 'T', 'N'};

    EXPECT_EQ(seqan3::to_char(dna5[0]), 'A');
    EXPECT_EQ(seqan3::to_char(dna5[1]), 'C');
    EXPECT_EQ(seqan3::to_char(dna5[2]), 'G');
    EXPECT_EQ(seqan3::to_char(dna5[3]), 'T');
    EXPECT_EQ(seqan3::to_char(dna5[4]), 'N');
}

TEST(dna5, assign_rank_to)
{
    EXPECT_EQ(seqan3::assign_rank_to(0, seqan::Dna5{}), seqan::Dna5{'A'});
    EXPECT_EQ(seqan3::assign_rank_to(1, seqan::Dna5{}), seqan::Dna5{'C'});
    EXPECT_EQ(seqan3::assign_rank_to(2, seqan::Dna5{}), seqan::Dna5{'G'});
    EXPECT_EQ(seqan3::assign_rank_to(3, seqan::Dna5{}), seqan::Dna5{'T'});
    EXPECT_EQ(seqan3::assign_rank_to(4, seqan::Dna5{}), seqan::Dna5{'N'});
}

TEST(dna5, assign_char_to)
{
    EXPECT_EQ(seqan3::assign_char_to('A', seqan::Dna5{}), seqan::Dna5{'A'});
    EXPECT_EQ(seqan3::assign_char_to('C', seqan::Dna5{}), seqan::Dna5{'C'});
    EXPECT_EQ(seqan3::assign_char_to('G', seqan::Dna5{}), seqan::Dna5{'G'});
    EXPECT_EQ(seqan3::assign_char_to('T', seqan::Dna5{}), seqan::Dna5{'T'});
    EXPECT_EQ(seqan3::assign_char_to('N', seqan::Dna5{}), seqan::Dna5{'N'});
}

TEST(dna5, complement)
{
    EXPECT_EQ(seqan3::complement(seqan::Dna5{'A'}), seqan::Dna5{'T'});
    EXPECT_EQ(seqan3::complement(seqan::Dna5{'C'}), seqan::Dna5{'G'});
    EXPECT_EQ(seqan3::complement(seqan::Dna5{'G'}), seqan::Dna5{'C'});
    EXPECT_EQ(seqan3::complement(seqan::Dna5{'T'}), seqan::Dna5{'A'});
    EXPECT_EQ(seqan3::complement(seqan::Dna5{'N'}), seqan::Dna5{'N'});
}

TEST(rna5, concepts)
{
    EXPECT_TRUE((seqan3::semialphabet<seqan::Rna5>));
    EXPECT_TRUE((seqan3::alphabet<seqan::Rna5>));
    EXPECT_TRUE((seqan3::nucleotide_alphabet<seqan::Rna5>));
}

TEST(rna5, to_rank)
{
    seqan::Rna5 rna5[] = {'A', 'C', 'G', 'T', 'N'};

    EXPECT_EQ(seqan3::to_rank(rna5[0]), 0);
    EXPECT_EQ(seqan3::to_rank(rna5[1]), 1);
    EXPECT_EQ(seqan3::to_rank(rna5[2]), 2);
    EXPECT_EQ(seqan3::to_rank(rna5[3]), 3);
    EXPECT_EQ(seqan3::to_rank(rna5[4]), 4);
}

TEST(rna5, to_char)
{
    seqan::Rna5 rna5[] = {'A', 'C', 'G', 'T', 'N'};

    EXPECT_EQ(seqan3::to_char(rna5[0]), 'A');
    EXPECT_EQ(seqan3::to_char(rna5[1]), 'C');
    EXPECT_EQ(seqan3::to_char(rna5[2]), 'G');
    EXPECT_EQ(seqan3::to_char(rna5[3]), 'U');
    EXPECT_EQ(seqan3::to_char(rna5[4]), 'N');
}

TEST(rna5, assign_rank_to)
{
    EXPECT_EQ(seqan3::assign_rank_to(0, seqan::Rna5{}), seqan::Rna5{'A'});
    EXPECT_EQ(seqan3::assign_rank_to(1, seqan::Rna5{}), seqan::Rna5{'C'});
    EXPECT_EQ(seqan3::assign_rank_to(2, seqan::Rna5{}), seqan::Rna5{'G'});
    EXPECT_EQ(seqan3::assign_rank_to(3, seqan::Rna5{}), seqan::Rna5{'T'});
    EXPECT_EQ(seqan3::assign_rank_to(4, seqan::Rna5{}), seqan::Rna5{'N'});
}

TEST(rna5, assign_char_to)
{
    EXPECT_EQ(seqan3::assign_char_to('A', seqan::Rna5{}), seqan::Rna5{'A'});
    EXPECT_EQ(seqan3::assign_char_to('C', seqan::Rna5{}), seqan::Rna5{'C'});
    EXPECT_EQ(seqan3::assign_char_to('G', seqan::Rna5{}), seqan::Rna5{'G'});
    EXPECT_EQ(seqan3::assign_char_to('T', seqan::Rna5{}), seqan::Rna5{'T'});
    EXPECT_EQ(seqan3::assign_char_to('N', seqan::Rna5{}), seqan::Rna5{'N'});
}

TEST(rna5, complement)
{
    EXPECT_EQ(seqan3::complement(seqan::Rna5{'A'}), seqan::Rna5{'T'});
    EXPECT_EQ(seqan3::complement(seqan::Rna5{'C'}), seqan::Rna5{'G'});
    EXPECT_EQ(seqan3::complement(seqan::Rna5{'G'}), seqan::Rna5{'C'});
    EXPECT_EQ(seqan3::complement(seqan::Rna5{'T'}), seqan::Rna5{'A'});
    EXPECT_EQ(seqan3::complement(seqan::Rna5{'N'}), seqan::Rna5{'N'});
}
