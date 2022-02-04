#pragma once

#include <seqan3/std/span>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace stellar::detail
{

template <typename sequence_t>
std::span<seqan3::dna4 const> reinterpret_as_dna4_sequence(sequence_t const & sequence)
{
    auto begin_it = std::ranges::begin(sequence);
    auto end_it = std::ranges::end(sequence);
    seqan3::dna4 const * begin_ptr = reinterpret_cast<seqan3::dna4 const *>(&*begin_it);
    seqan3::dna4 const * end_ptr = reinterpret_cast<seqan3::dna4 const *>(&*end_it);

    return std::span<seqan3::dna4 const>{begin_ptr, end_ptr};
}

} // namespace stellar::detail
