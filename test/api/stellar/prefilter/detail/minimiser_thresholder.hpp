
#pragma once

#include <stellar/prefilter/detail/valik_local_query_prefilter_index.hpp>

struct minimiser_thresholder
{
    minimiser_thresholder(
        stellar::detail::valik_local_query_prefilter_index const & lqpindex,
        size_t const minLength,
        size_t const maxError)
    {
        auto search_arguments = stellar::detail::valik_local_prefilter_search_arguments(lqpindex, minLength, maxError);

        threshold = stellar::detail::local_prefilter_threshold{search_arguments};
    }

    template <typename count_t>
    std::vector<bool> bin_thresholds_met(seqan3::counting_vector<count_t> const & bin_minimiser_counts, size_t const minimiser_count)
    {
        std::vector<bool> thresholds{};

        size_t const minimiser_threshold = get(minimiser_count);

        for (size_t bin_minimiser_count : bin_minimiser_counts)
        {
            thresholds.push_back(bin_minimiser_count >= minimiser_threshold);
        }

        return thresholds;
    }

    size_t get(size_t const minimiser_count) const
    {
        return threshold.get(minimiser_count);
    }

    size_t min_number_of_minimisers() const
    {
        return threshold.min_number_of_minimisers;
    }

    size_t max_number_of_minimisers() const
    {
        return threshold.max_number_of_minimisers;
    }

    stellar::detail::local_prefilter_threshold threshold;
};
