
#pragma once

#include <stellar/test/alignment_fragment.hpp>

namespace stellar::test::fixture::a512_simSeq_0dot05::forward
{
// database file: 512_simSeq1_5e-2.fa
// query file: 512_simSeq2_5e-2.fa
// minimal match length             : 150
// maximal error rate (epsilon)     : 0.05
// maximal x-drop                   : 10
// k-mer (q-gram) length            : 7
// search forward strand            : yes
// search reverse complement        : no

// best match
// seq1	Stellar	eps-matches	203842	204051	95.3051	+	.	seq2;seq2Range=797542,797749;cigar=41M1D48M1D39M1I49M1D5M1I17M1D1M1D3M1I2M;mutations=27A,129A,148T,184C,206A
static stellar::test::alignment_fragment<seqan::Dna5> const fragment0
{
    "TTTGAAAAAACTTACTTGTCAATTT[",
    "AACGGGGTTGAACCAATTATTCATC[",

    // 0 ... 100
    " ATAAGAAGCGGGAATTTTGTATTGCTTCAATACAGTCTGTCGTTTGGGCACGACCAGTTCGGACCTCGCACCAATCGCTGACGATGCGGCAACGTACACA",
    //|||||||||||||||||||||||||| |||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||| |||||||||
    " ATAAGAAGCGGGAATTTTGTATTGCTACAATACAGTCTGTC-TTTGGGCACGACCAGTTCGGACCTCGCACCAATCGCTGACGATGCGGC-ACGTACACA",

    // 100 ... 200
    " AACATGCGCATCTAACCGCTAGCTCACTAT-GACTATTTCGCTGATTATCTGGTGGGGGCATCCACTATCCGTATTGGCCGAAATA-CAGCGGCGTGATG",
    //|||||||||||||||||||||||||||||| |||||||||||||||||| |||||||||||||||||||||||||||||| ||||| |||||||||||||
    " AACATGCGCATCTAACCGCTAGCTCACTATAGACTATTTCGCTGATTATTTGGTGGGGGCATCCACTATCCGTATTGGCC-AAATACCAGCGGCGTGATG",

    // 200 ... 213
    " AATCACGTCG-GG",
    //|||| | ||| ||
    " AATC-C-TCGAGG",

    "]AGTAGACGCCTCTCTTGTACTCGGG",
    "]TTGTGGCATAAGCGAGCATACAGCG"
};

} // namespace stellar::test::fixture::a512_simSeq_0dot05::forward
