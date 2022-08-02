#pragma once

#include <algorithm>
#include <cstddef>
#include <deque>

#include <stellar/test/prefilter/local_query_prefilter/minimiser_window/minimiser_state.hpp>

namespace stellar {
namespace test {

template <typename value_type>
struct seqan3_minimiser_window
{
    seqan3_minimiser_window(size_t const window_size)
        : window_size{window_size}
    {}

    value_type min() const
    {
        return minimiser_value;
    }

    //!\brief Calculates minimisers for the first window.
    template <typename iterator_t, typename sentinel_t>
    iterator_t initialize(iterator_t iterator, sentinel_t sentinel)
    {
        size_t const size = sentinel - iterator;
        size_t const window_size = std::min<size_t>(this->window_size, size);

        for (size_t i = 0u; i < window_size - 1u; ++i)
        {
            window_values.push_back(*iterator);
            ++iterator;
        }
        window_values.push_back(*iterator);

        recalculate_minimum();

        return iterator;
    }

    /*!\brief Calculates the next minimiser value.
     * \returns True, if new minimiser is found or end is reached. Otherwise returns false.
     * \details
     * For the following windows, we remove the first window value (is now not in window_values) and add the new
     * value that results from the window shifting.
     */
    minimiser_state cyclic_push(value_type const new_value)
    {
        window_values.pop_front();
        window_values.push_back(new_value);

        if (minimiser_position_offset == 0)
        {
            recalculate_minimum();
            return minimiser_state::left_window;
        }

        if (new_value < minimiser_value)
        {
            minimiser_value = new_value;
            minimiser_position_offset = window_values.size() - 1;
            return minimiser_state::new_minimizer;
        }

        --minimiser_position_offset;
        return minimiser_state::unchanged;
    }

protected:

    void recalculate_minimum()
    {
        auto minimiser_it = std::min_element(window_values.begin(), window_values.end(), std::less_equal<value_type>{});
        minimiser_value = *minimiser_it;
        minimiser_position_offset = std::distance(std::begin(window_values), minimiser_it);
    }

    //!\brief The minimiser value.
    value_type minimiser_value{};

    //!\brief The offset relative to the beginning of the window where the minimizer value is found.
    size_t minimiser_position_offset{};

    //!\brief Stored values per window. It is necessary to store them, because a shift can remove the current minimiser.
    std::deque<value_type> window_values{};

    size_t window_size;
};

} // namespace stellar::test
} // namespace stellar
