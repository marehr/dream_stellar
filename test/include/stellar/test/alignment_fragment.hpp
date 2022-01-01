
#pragma once

#include <string_view>
#include <seqan3/std/span>

#include <stellar/stellar.hpp>

namespace stellar
{
namespace test
{

template <typename TAlphabet>
struct alignment_fragment
{
    alignment_fragment(alignment_fragment const &) = default/*delete*/;
    alignment_fragment & operator=(alignment_fragment const &) = default/*delete*/;

    alignment_fragment(std::initializer_list<std::string_view const> aligned_sequences)
    {
        resize(rows(alignment), 2);

        if ((aligned_sequences.size() % 2) != 0u)
        {
            throw std::runtime_error{"To few arguments."};
        }

        std::vector<std::string_view> aligned_sequences1{};
        std::vector<std::string_view> aligned_sequences2{};

        bool even = true;
        for (std::string_view const & aligned_sequence: aligned_sequences)
        {
            std::vector<std::string_view> & active = even ? aligned_sequences1 : aligned_sequences2;
            active.push_back(aligned_sequence);
            even = !even;
        }

        parse_aligned_sequence(std::span<std::string_view const>{aligned_sequences1},
                               sequence1, beginPosition1, endPosition1, row(alignment, 0));
        parse_aligned_sequence(std::span<std::string_view const>{aligned_sequences2},
                               sequence2, beginPosition2, endPosition2, row(alignment, 1));
    }

    using TSequence = seqan::String<TAlphabet>;
    using TAlign = seqan::Align<TSequence const, seqan::ArrayGaps>;

    TSequence sequence1{};
    TSequence sequence2{};

    size_t beginPosition1{};
    size_t beginPosition2{};

    size_t endPosition1{};
    size_t endPosition2{};

    TAlign alignment{};

private:
    template <typename TGappedSequence>
    void parse_aligned_sequence(
        std::span<std::string_view const> const splitAlignedSequenceWithWhitespace,
        TSequence & sequence,
        size_t & beginPosition,
        size_t & endPosition,
        TGappedSequence & gappedSequence)
    {
        auto appendPureSequenceTo = [](std::string_view const segment, auto & sequence, auto && filterOut)
        {
            for (char const chr : segment)
            {
                if (filterOut(chr))
                    continue;

                appendValue(sequence, chr);
            }
        };

        // concatenate all sequences together and remove whitespaces
        std::vector<char> alignedSequenceVector;
        for (std::string_view const & alignedSequenceWithWhitespace: splitAlignedSequenceWithWhitespace)
        {
            appendPureSequenceTo(alignedSequenceWithWhitespace, alignedSequenceVector, [](char const chr)
            {
                return chr == ' ';
            });
        }
        std::string_view const alignedSequence{alignedSequenceVector.data(), alignedSequenceVector.size()};

        // copy sequence without gap symbols or start/end marker
        appendPureSequenceTo(alignedSequence, sequence, [](char const chr)
        {
            return chr == '-' || chr == '[' || chr == ']';
        });

        setSource(gappedSequence, sequence);

        // find start position of '['
        beginPosition = alignedSequence.find('[');
        assert(beginPosition != alignedSequence.npos); // sequence must include "[" char

        // find end position of ']' in aligned sequence
        size_t const alignedEndPosition = alignedSequence.find(']', beginPosition);
        assert(alignedEndPosition != alignedSequence.npos); // sequence must include "]" char

        // count gaps in between '[', ... ,']';
        // aligned sequence length minus 2 markers minus unaligned sequence length
        size_t const gapCount = alignedSequence.size() - 2u - length(sequence);

        // over counted alignment start marker '[' and gap symbols '-'.
        endPosition = alignedEndPosition - gapCount - 1;

        // set begin and end position
        setBeginPosition(gappedSequence, beginPosition);
        setEndPosition(gappedSequence, endPosition);

        // copy over alignment
        size_t i = 0;
        for (char const chr : alignedSequence.substr(beginPosition + 1, endPosition - beginPosition + gapCount))
        {
            bool const is_gap = chr == '-';
            if (is_gap)
            {
                insertGap(gappedSequence, i);
            } else
            {
                // Note: we first need to convert to TAlphabet, otherwise
                // we get some SegFault due to some initialization order
                // shenanigans if we use this function in a global initializer.
                EXPECT_EQ(chr, convert<TAlphabet>(gappedSequence[i])) << "[i]: " << i;
            }
            EXPECT_EQ(isCharacter(gappedSequence, i), !is_gap) << "[i]: " << i;
            EXPECT_EQ(isGap(gappedSequence, i), is_gap) << "[i]: " << i;

            ++i;
        }
    }
};

} // namespace stellar::test
} // namespace stellar
