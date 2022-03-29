#pragma once

#include <valik/shared.hpp> // search_arguments

#include <stellar/prefilter/detail/valik_local_query_prefilter_index.hpp>

namespace stellar::detail
{

inline valik::search_arguments
valik_local_prefilter_search_arguments(
    stellar::detail::valik_local_query_prefilter_index const & lqpindex,
    size_t const minLength,
    size_t const maxError)
{
    seqan3::shape ungapped_shape{seqan3::ungapped{lqpindex._kmer_size}};
    valik::search_arguments const arguments
    {
        .window_size{lqpindex._window_size.get()},
        .shape{ungapped_shape}, // kmer_size
        .shape_size{ungapped_shape.size()}, // TODO: really needed?
        .shape_weight{ungapped_shape.count()}, // TODO: really needed?
        // step size is i + pattern_size - overlap = i + 1
        // => pattern_size - overlap = 1
        // => overlap = pattern_size - 1
        // TODO: checking every possible overlap is to much and there must be a smarter way
        .overlap{minLength - 1u},
        .fpr{0.00},
        .errors{static_cast<uint8_t>(maxError)},
        // this is the max-size how many minimiser need to match a threshold?!
        // it splits the query sequence into many chunks of size pattern_size
        .pattern_size{minLength} // >= window_size
    };

    return arguments;
}

} // namespace stellar::detail
