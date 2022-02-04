#pragma once

#include <valik/shared.hpp> // search_arguments
#include <valik/search/compute_simple_model.hpp> // threshold

#ifdef STELLAR_MINIMISER_THRESHOLD_INLINE_SRC
// include_dir = /home/marehr/.../dream_stellar/lib/valik/include/
// relative from include_dir (valik/include)
//     ../src/search/minimiser_model.cpp
//     ../src/search/detail/helper.cpp
#include "../src/search/minimiser_model.cpp"
#include "../src/search/detail/helper.cpp"
#endif // STELLAR_MINIMISER_THRESHOLD_INLINE_SRC

namespace stellar::detail
{

struct local_prefilter_threshold : public valik::threshold
{
    local_prefilter_threshold() = default;

    local_prefilter_threshold(valik::search_arguments const & search_arguments) :
        valik::threshold{valik::make_threshold_data(search_arguments)},
        treshold_was_set{search_arguments.treshold_was_set},
        threshold{search_arguments.threshold}
    {
        // current local_prefilter uses +2 and doesn't count one minimiser
        for (size_t & threshold: this->precomp_thresholds)
            threshold = std::max<size_t>(threshold, 3u) - 3;
    }

    size_t get(size_t const minimiser_count) const
    {
        auto const & threshold_data = *this;
        auto const & arguments = *this;
        return
            arguments.treshold_was_set ?
            static_cast<size_t>(minimiser_count * arguments.threshold) : threshold_data.kmers_per_window == 1 ?
            threshold_data.kmer_lemma : threshold_data.precomp_thresholds[std::min(minimiser_count < threshold_data.min_number_of_minimisers ?
            0 : minimiser_count - threshold_data.min_number_of_minimisers,
            threshold_data.max_number_of_minimisers - threshold_data.min_number_of_minimisers)] + 2;
    }

    bool treshold_was_set{false};
    double threshold{};
};

inline local_prefilter_threshold
valik_local_prefilter_make_threshold(
    valik::search_arguments const & search_arguments)
{
    // make_threshold_data
    return local_prefilter_threshold{search_arguments};
}

} // namespace stellar::detail
