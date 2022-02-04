#pragma once

#include <seqan3/alphabet/range/sequence.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>

namespace stellar
{

template <typename sequence_t>
concept reverse_complementable_sequence =
    seqan3::sequence<sequence_t> && // sequence
    seqan3::nucleotide_alphabet<std::ranges::range_reference_t<sequence_t>> && // complementable
    std::ranges::bidirectional_range<sequence_t>; // reversable
} // namespace stellar
