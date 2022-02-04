#pragma once

// DEBUG: remove me
#include <seqan3/alphabet/views/to_rank.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <bitset>

#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <valik/search/local_prefilter.hpp>
#include <valik/search/search_setup.hpp>

#ifdef STELLAR_MINIMISER_THRESHOLD_INLINE_SRC
// include_dir = /home/marehr/.../dream_stellar/lib/sliding-window/include/
// relative from include_dir (valik/include)
//     ../src/search/minimiser_model.cpp
//     ../src/search/detail/helper.cpp
#include "../src/search/minimiser_model.cpp"
#include "../src/search/detail/helper.cpp"
#endif // STELLAR_MINIMISER_THRESHOLD_INLINE_SRC

#include <stellar/prefilter/detail/reinterpret_as_dna4_sequence.hpp>
#include <stellar/prefilter/detail/valik_local_query_prefilter_index.hpp>

namespace stellar::detail
{

template <typename TAlphabet, typename TOnBinResultFn>
void valik_local_prefilter(
    seqan::String<TAlphabet> const & sequence,
    stellar::detail::valik_local_query_prefilter_index const & lqpindex,
    size_t minLength,
    TOnBinResultFn && on_bin_result)
{
    valik_local_prefilter(sequence, lqpindex, minLength, 0u, on_bin_result);
}

template <typename TAlphabet, typename TOnBinResultFn>
void valik_local_prefilter(
    seqan::String<TAlphabet> const & sequence,
    stellar::detail::valik_local_query_prefilter_index const & lqpindex,
    size_t minLength,
    size_t maxError,
    TOnBinResultFn && on_bin_result)
{
    assert(minLength >= 1u);

    if (seqan::length(sequence) < minLength)
    {
        std::cerr <<
            "WARNING: sequence length (= " << seqan::length(sequence) << ") "
            "is smaller than the minLength (= " << minLength << "). "
            "Skipping this sequence." << std::endl;
        return;
    }

    if (minLength < lqpindex._window_size.get())
    {
        std::cerr <<
            "WARNING: minLength (= " << minLength << ") "
            "is smaller than the windowSize (= " << lqpindex._window_size.get() << ") "
            "of the minimiser. "
            "The prefilter can't be used and thus this query will be searched in all databases." << std::endl;

        for (size_t bin_id = 0u, bin_count = lqpindex._ibf.bin_count(); bin_id < bin_count; ++bin_id)
        {
            on_bin_result(bin_id);
        }

        return;
    }

    std::span<seqan3::dna4 const> dna4_sequence_view = stellar::detail::reinterpret_as_dna4_sequence(sequence);

    valik::query_record record{};
    record.sequence = std::vector<seqan3::dna4>{dna4_sequence_view.begin(), dna4_sequence_view.end()};

#if STELLAR_DEBUG_MINIMISER
    std::cout << "sequence: " << sequence << std::endl;
    seqan3::debug_stream << "rank_sequence: " << (sequence | seqan3::views::to_rank) << std::endl;
    seqan3::debug_stream << "dna4_sequence: " << record.sequence << std::endl;
    seqan3::debug_stream << "dna4_rank_sequence: " << (record.sequence | seqan3::views::to_rank) << std::endl;

    std::cout << "lqpindex._kmer_size: " << (int)lqpindex._kmer_size.value << std::endl;
    std::cout << "lqpindex._window_size: " << (int)lqpindex._window_size.get() << std::endl;
    std::cout << "lqpindex._seed: " << (int)lqpindex._seed.get() << std::endl;

    print_all_minimiser_with_containment(sequence, lqpindex);
#endif // STELLAR_DEBUG_MINIMISER

    std::span<valik::query_record const> const records{&record, 1u};
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> const & ibf = lqpindex._ibf;

    seqan3::shape ungapped_shape{seqan3::ungapped{lqpindex._kmer_size}};
    valik::search_arguments const arguments
    {
        .window_size{lqpindex._window_size.get()},
        .shape{ungapped_shape}, // kmer_size
        .shape_size{ungapped_shape.size()}, // TODO: really needed?
        .shape_weight{ungapped_shape.count()}, // TODO: really needed?
        // this is the max-size how many minimiser need to match a threshold?!
        // it splits the query sequence into many chunks of size pattern_size
        .pattern_size{minLength}, // >= window_size
        // step size is i + pattern_size - overlap = i + 1
        // => pattern_size - overlap = 1
        // => overlap = pattern_size - 1
        // TODO: checking every possible overlap is to much and there must be a smarter way
        .overlap{minLength - 1u},
        .errors{static_cast<uint8_t>(maxError)}
    };

    // make_threshold_data
    valik::threshold threshold_data = valik::make_threshold_data(arguments);
    // current local_prefilter uses +2 and doesn't count one minimiser
    for (size_t & threshold: threshold_data.precomp_thresholds)
        threshold = std::max<size_t>(threshold, 3u) - 3;

    std::vector<valik::query_result> query_results = valik::local_prefilter(
        /*std::span<query_record const> const &*/ records,
        /*seqan3::interleaved_bloom_filter<ibf_data_layout> const &*/ ibf,
        /*search_arguments const &*/ arguments,
        /*threshold const &*/ threshold_data
    );

    assert(query_results.size() > 0u);

    for (valik::query_result const & query_result : query_results)
    {
        for (size_t bin_id : query_result.get_hits())
            on_bin_result(bin_id);
    }
}
} // namespace stellar::detail
