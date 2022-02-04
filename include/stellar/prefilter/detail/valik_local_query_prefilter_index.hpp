#pragma once

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <valik/shared.hpp>

#include <stellar/seqan2to3_alphabet.hpp>
#include <stellar/prefilter/detail/reinterpret_as_dna4_sequence.hpp>
#include <stellar/prefilter/detail/reverse_complementable_sequence.hpp>

namespace stellar::detail
{

struct valik_local_query_prefilter_index
{
    static constexpr seqan3::data_layout data_layout_mode = seqan3::data_layout::uncompressed;

    inline constexpr static uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
    {
        return valik::adjust_seed(kmer_size, seed);
    }

    template <typename sequence_t>
        requires (stellar::reverse_complementable_sequence<sequence_t>)
    auto hash_view(sequence_t & sequence) const
    {
        return seqan3::views::minimiser_hash(
            reinterpret_as_dna4_sequence(sequence),
            seqan3::ungapped{_kmer_size},
            seqan3::window_size{_window_size},
            seqan3::seed{_seed});
    }

    template <std::ranges::viewable_range sequence_t>
        requires (stellar::reverse_complementable_sequence<sequence_t>)
    auto insert_sequence(seqan3::bin_index const bin_index, sequence_t && sequence)
    {
        auto minimiser_view = hash_view(sequence);
        auto const minimiser_begin = minimiser_view.begin();
        auto const minimiser_end = minimiser_view.end();
        auto minimiser_it = minimiser_begin;

        for (; minimiser_it != minimiser_end; ++minimiser_it)
        {
            size_t const hash = *minimiser_it;
            _ibf.emplace(hash, bin_index);
        }
    }

    seqan3::interleaved_bloom_filter<data_layout_mode> _ibf{};

    seqan3::ungapped _kmer_size{};
    seqan3::window_size _window_size{};
    seqan3::seed _seed{adjust_seed(_kmer_size.value)};
};


} // namespace stellar::detail
