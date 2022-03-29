#pragma once

#include <valik/shared.hpp> // search_arguments
#include <raptor/threshold/threshold.hpp> // threshold

namespace stellar::detail
{

struct local_prefilter_threshold : public raptor::threshold::threshold
{
    using base_t = raptor::threshold::threshold;

    local_prefilter_threshold() = default;

    local_prefilter_threshold(valik::search_arguments const & search_arguments) :
        base_t{search_arguments.make_threshold_parameters()}//,
        // treshold_was_set{search_arguments.treshold_was_set},
        // threshold{search_arguments.threshold}
    {
        // // current local_prefilter uses +2 and doesn't count one minimiser
        // for (size_t & threshold: this->precomp_thresholds)
        //     threshold = std::max<size_t>(threshold, 3u) - 3;
    }

    size_t get(size_t const minimiser_count) const
    {
        base_t const & threshold = *this;
        return threshold.get(minimiser_count);
    }

    // bool treshold_was_set{false};
    // double threshold{};
};

inline local_prefilter_threshold
valik_local_prefilter_make_threshold(
    valik::search_arguments const & search_arguments)
{
    // make_threshold_data
    return local_prefilter_threshold{search_arguments};
}

} // namespace stellar::detail
