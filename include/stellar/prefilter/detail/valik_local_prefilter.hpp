#pragma once

#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <valik/search/local_prefilter.hpp>
#include <valik/search/search_setup.hpp>

#include <stellar/prefilter/detail/reinterpret_as_dna4_sequence.hpp>
#include <stellar/prefilter/detail/valik_local_prefilter_search_arguments.hpp>
#include <stellar/prefilter/detail/valik_local_prefilter_threshold.hpp>
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

    std::span<valik::query_record const> const records{&record, 1u};
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> const & ibf = lqpindex._ibf;

    auto search_arguments = valik_local_prefilter_search_arguments(lqpindex, minLength, maxError);
    auto threshold = valik_local_prefilter_make_threshold(search_arguments);

    std::vector<valik::query_result> query_results = valik::local_prefilter(
        /*std::span<query_record const> const &*/ records,
        /*seqan3::interleaved_bloom_filter<ibf_data_layout> const &*/ ibf,
        /*search_arguments const &*/ search_arguments,
        /*threshold const &*/ threshold
    );

    assert(query_results.size() > 0u);

    for (valik::query_result const & query_result : query_results)
    {
        for (size_t bin_id : query_result.get_hits())
            on_bin_result(bin_id);
    }
}
} // namespace stellar::detail
