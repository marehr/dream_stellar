#pragma once

#if STELLAR_DEBUG_MINIMISER
// TODO: debug remove me
#include <bitset>
#include <seqan3/core/debug_stream.hpp>

#include <stellar/prefilter/detail/reinterpret_as_dna4_sequence.hpp>

namespace stellar::detail
{

template <typename minimiser_it_t, typename underlying_sequence_t>
void print_minimiser(
    minimiser_it_t const & minimiser_it,
    minimiser_it_t const & minimiser_begin_it,
    size_t const kmer_size,
    underlying_sequence_t const & underlying_sequence,
    std::string_view sequence_name = "sequence")
{
    size_t const hash = *minimiser_it;

    // in hash coordinates
    size_t const hash_window_begin_position = minimiser_it.base() - minimiser_begin_it.base();
    size_t const hash_window_size = minimiser_it.window_size();
    size_t const hash_minimiser_offset = minimiser_it.minimiser_offset();

    // in sequence coordinates
    size_t const window_size = hash_window_size + kmer_size - 1;
    size_t const window_begin_position = hash_window_begin_position;
    size_t const minimiser_begin_position = window_begin_position + hash_minimiser_offset;
    size_t const minimiser_end_position = minimiser_begin_position + kmer_size;
    size_t const window_end_position = window_begin_position + window_size;

    auto infix = [&](size_t const begin_position, size_t const end_position)
    {
        std::span<seqan3::dna4 const> span = reinterpret_as_dna4_sequence(underlying_sequence);

        return span.subspan(begin_position, end_position - begin_position);
    };

    seqan3::debug_stream
        << sequence_name << ": "
        << infix(0u, window_begin_position)
        << "["
        << infix(window_begin_position, minimiser_begin_position)
        << "|"
        << infix(minimiser_begin_position, minimiser_end_position)
        << "|"
        << infix(minimiser_end_position, window_end_position)
        << "]"
        << infix(window_end_position, std::ranges::size(underlying_sequence))
        << " hash(" << hash << ")";
}

template <typename sequence_t, typename hash_adaptor_fn>
void print_all_minimiser(sequence_t && sequence, size_t const kmer_size, hash_adaptor_fn && hash_adaptor)
{
    auto minimiser_view = hash_adaptor(sequence);
    auto const minimiser_begin = minimiser_view.begin();
    auto const minimiser_end = minimiser_view.end();
    auto minimiser_it = minimiser_begin;
    for (; minimiser_it != minimiser_end; ++minimiser_it)
    {
        print_minimiser(minimiser_it, minimiser_begin, kmer_size, sequence);
        std::cout << std::endl;
    }
}

template <typename sequence_t, typename lqpindex_t>
void print_all_minimiser_with_containment(
    sequence_t && sequence,
    lqpindex_t const & lqpindex,
    std::string_view sequence_name = "query[]")
{
    auto agent = lqpindex._ibf.membership_agent();

    auto minimiser_view = lqpindex.hash_view(sequence);
    auto const minimiser_begin = minimiser_view.begin();
    auto const minimiser_end = minimiser_view.end();
    auto minimiser_it = minimiser_begin;
    for (size_t i = 0; minimiser_it != minimiser_end; ++minimiser_it, ++i)
    {
        size_t const minimiser = *minimiser_it;
        size_t bit_sequence = agent.bulk_contains(minimiser).raw_data().get_int(0u);

        seqan3::debug_stream << i << ". ";
        print_minimiser(minimiser_it, minimiser_begin, lqpindex._kmer_size.value, sequence, sequence_name);
        seqan3::debug_stream << " contained in "
                             << std::bitset<5u>{bit_sequence};
        seqan3::debug_stream << std::endl;
    }
}

} // namespace stellar::detail
#endif // STELLAR_DEBUG_MINIMISER
