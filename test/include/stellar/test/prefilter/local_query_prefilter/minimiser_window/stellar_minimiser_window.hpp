#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

namespace stellar {
namespace test {

template <typename value_type>
struct stellar_minimiser_window
{
    // Memory Layout
    // s1, s2, ..., sW: sorted elements
    // e: one element between sorted and unsorted
    // u1, u2, ..., uW: unsorted elements
    // [ s1, s2, ..., sW, e, u1, u2, ..., uW]
    stellar_minimiser_window(size_t const window_size) :
        window_size{window_size}
    {
        window_values.resize(2 * window_size + 1);
    }

    value_type const & min() const
    {
        return *minimiser_it;
    }

    //!\brief Calculates minimisers for the first window.
    template <typename iterator_t, typename sentinel_t>
    iterator_t initialize(iterator_t input_iterator, sentinel_t input_sentinel)
    {
        size_t const size = input_sentinel - input_iterator;
        size_t const window_size = std::min<size_t>(this->window_size, size);
        assert(window_size > 0u);

        this->sorted_begin = window_values.data();
        this->sorted_end = window_values.data() + this->window_size;

        iterator_t new_input_sentinel = input_iterator + window_size;
        this->unsorted_begin = this->sorted_end + 1;
        this->unsorted_end = std::copy(input_iterator, new_input_sentinel, this->unsorted_begin);

        // initial minimiser is the last element in the sorted range
        this->minimiser_it = this->sorted_end - 1;

        recalculate_minimum();

        return --new_input_sentinel;
    }

    /*!\brief Calculates the next minimiser value.
     * \returns True, if new minimiser is found or end is reached. Otherwise returns false.
     * \details
     * For the following windows, we remove the first window value (is now not in window_values) and add the new
     * value that results from the window shifting.
     */
    bool cyclic_push(value_type const new_value)
    {
        this->diagnostics();
        ++this->sorted_begin;

        *this->unsorted_end = new_value;
        this->unsorted_minimiser_it = indexed_minimum(this->unsorted_minimiser_it, this->unsorted_end);
        ++this->unsorted_end;

        value_type * previous_minimiser_it = this->minimiser_it;
        if (this->sorted_begin == this->sorted_end)
        {
            std::cout << "REBUILD!" << std::endl;
            std::cout << "this->minimiser_it: " << (this->minimiser_it - window_values.data()) << std::endl;
            std::cout << "this->unsorted_minimiser_it: " << (this->unsorted_minimiser_it - window_values.data()) << std::endl;
            // TODO! not every rebuild changes minimum
            // initial minimiser is either the last element in the sorted range,
            // or if the current minimiser is in the unsorted list it must stay the same
            bool minimiser_in_unsorted = this->minimiser_it == this->unsorted_minimiser_it;
            std::cout << "minimiser_in_unsorted: " << (minimiser_in_unsorted ? "TRUE" : "FALSE") << std::endl;
            previous_minimiser_it = minimiser_in_unsorted ? this->unsorted_minimiser_it - this->window_size : nullptr;
            // Note: that in recalculate_minimum copies all unsorted elements into the sorted area
            // The sorted area is offset by this->window_size from the unsorted area.
            this->minimiser_it = minimiser_in_unsorted ? previous_minimiser_it : this->sorted_end - 1;
            recalculate_minimum();
        } else
        {
            bool const minimiser_left_window = this->minimiser_it < this->sorted_begin;
            this->minimiser_it = (minimiser_left_window ? this->sorted_begin : this->minimiser_it);
            this->minimiser_it = indexed_minimum(this->minimiser_it, this->unsorted_minimiser_it);
        }

        return previous_minimiser_it != this->minimiser_it;
    }

protected:

    void diagnostics()
    {
        std::cout << "stellar_minimiser_window:" << std::endl;
        std::cout << "\tsorted:   [";
        for (auto it = this->sorted_begin; it != this->sorted_end; ++it)
            std::cout << (it == this->minimiser_it ? "*" : "") << *it << ", ";
        std::cout << "]" << std::endl;
        std::cout << "\tunsorted: [";
        for (auto it = this->unsorted_begin; it != this->unsorted_end; ++it)
            std::cout << (it == this->minimiser_it ? "*" : "") << (it == this->unsorted_minimiser_it ? "!" : "") << *it << ", ";
        std::cout << "]" << std::endl;
    }

    static constexpr value_type * indexed_minimum(value_type * const old_it, value_type * const new_it)
    {
        return *new_it < *old_it ? new_it : old_it;
    }

    void recalculate_minimum()
    {
        assert(this->sorted_end == this->unsorted_begin);

        value_type * sorted_it = this->sorted_end;
        value_type * unsorted_it = this->unsorted_end;
        value_type * unsorted_sentinel = this->unsorted_begin;

        // unsorted is non-empty
        assert(unsorted_it != unsorted_sentinel);

        // update pointer to be e in memory layout
        // [ s1, s2, ..., sW, e, u1, u2, ..., uW]
        this->unsorted_minimiser_it = sorted_it;

        // last element is minimum
        // NOTE: it is also the minimum for the future unsorted values
        // i.e. if never a smaller value comes, it will always be the "global" minimum
        --unsorted_it;
        --sorted_it;
        *sorted_it = *unsorted_it;

        // we need to update value to be the value of the last element
        *this->unsorted_minimiser_it = *sorted_it;

        // construct minimum from last to first element
        for (; unsorted_it != unsorted_sentinel; )
        {
            --unsorted_it;
            --sorted_it;
            *sorted_it = *unsorted_it;
            this->minimiser_it = indexed_minimum(this->minimiser_it, sorted_it);
            *sorted_it = *this->minimiser_it;
        }

        // update pointer
        this->sorted_begin = sorted_it; // sorted region is now non-empty
        this->unsorted_end = this->unsorted_begin; // unsorted region is now empty
    }

    value_type * minimiser_it{};
    value_type * unsorted_minimiser_it{};

    //!\brief Stored values per window. It is necessary to store them, because a shift can remove the current minimiser.
    std::vector<value_type> window_values{};

    value_type * sorted_begin;
    value_type * sorted_end;

    value_type * unsorted_begin;
    value_type * unsorted_end;

    size_t window_size;
};

} // namespace stellar::test
} // namespace stellar
