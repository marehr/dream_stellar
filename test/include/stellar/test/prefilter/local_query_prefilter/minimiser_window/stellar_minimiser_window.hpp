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
    // e*: one element between sorted and unsorted
    // u1, u2, ..., uW: unsorted elements
    // [ s1, s2, ..., sW, e*, u1, u2, ..., uW]
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

        this->sorted_begin = {window_values.data(), this};
        this->sorted_end = this->sorted_begin + this->window_size;

        iterator_t new_input_sentinel = input_iterator + window_size;
        this->unsorted_begin = {this->sorted_end.ptr + 1, this};
        this->unsorted_end = {std::copy(input_iterator, new_input_sentinel, this->unsorted_begin.ptr), this};

        // initial minimiser is the first element in the sorted range
        //
        // Note recalculate_minimum will do:
        // In range (sorted_end, minimiser_it] recompute minimiser backwards (old minimiser would have left window)
        // In range (minimiser_it, sorted_begin] set old minimiser value as it is still active in that range
        //
        // Setting minimiser_it = sorted_begin makes the second range empty.
        this->minimiser_it = (mixed_ptr)(this->sorted_end - 1);
        this->unsorted_minimiser_it = this->unsorted_begin;

        std::cout << "INIT!" << std::endl;
        recalculate_minimum();

        this->diagnostics();

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
        ++this->sorted_begin;

        mixed_ptr previous_minimiser_it = this->minimiser_it;

        { // append new element
            *this->unsorted_end = new_value;
            this->unsorted_minimiser_it = indexed_minimum(this->unsorted_minimiser_it, this->unsorted_end);
            ++this->unsorted_end;

            bool const minimiser_left_window = this->minimiser_it < (mixed_ptr)this->sorted_begin;
            this->minimiser_it = (minimiser_left_window ? (mixed_ptr)(this->sorted_begin) : this->minimiser_it);
            // TODO: not sure if correct, but we filter out same value runs
            bool is_same_sorted_value = minimiser_left_window && *previous_minimiser_it == *this->minimiser_it;
            previous_minimiser_it = is_same_sorted_value ? this->minimiser_it : previous_minimiser_it;
            // this->unsorted_minimiser_it = this->unsorted_minimiser_it;
            this->minimiser_it = indexed_minimum<mixed_ptr>(this->minimiser_it, (mixed_ptr)this->unsorted_minimiser_it);
        }
        this->diagnostics();

        bool minimiser_changed = previous_minimiser_it != this->minimiser_it;

        if (this->sorted_begin == this->sorted_end)
        {
            std::cout << "REBUILD!" << std::endl;
            std::cout << "previous_minimiser_it: " << (previous_minimiser_it - mixed_ptr{window_values.data(), this}) << std::endl;
            std::cout << "this->minimiser_it: " << (this->minimiser_it - mixed_ptr{window_values.data(), this}) << std::endl;
            std::cout << "this->unsorted_minimiser_it: " << (mixed_ptr{this->unsorted_minimiser_it} - mixed_ptr{window_values.data(), this}) << std::endl;
            // TODO! not every rebuild changes minimum
            // initial minimiser is either the last element in the sorted range,
            // or if the current minimiser is in the unsorted list it must stay the same
            // bool minimiser_in_unsorted = this->minimiser_it == (mixed_ptr)this->unsorted_minimiser_it;
            // std::cout << "minimiser_in_unsorted: " << (minimiser_in_unsorted ? "TRUE" : "FALSE") << std::endl;
            // // Note: that in recalculate_minimum copies all unsorted elements into the sorted area
            // // The sorted area is offset by this->window_size from the unsorted area.
            // this->minimiser_it =
            //     minimiser_in_unsorted ?
            //     mixed_ptr{this->unsorted_minimiser_it} - (this->window_size + 1) :
            //     mixed_ptr{this->sorted_begin};
            // std::cout << "NEW this->minimiser_it: " << (this->minimiser_it - mixed_ptr{window_values.data(), this}) << std::endl;
            recalculate_minimum();

            this->diagnostics();
            std::cout << "AFTER: previous_minimiser_it: " << (previous_minimiser_it - mixed_ptr{window_values.data(), this}) << std::endl;
            std::cout << "AFTER: this->minimiser_it: " << (this->minimiser_it - mixed_ptr{window_values.data(), this}) << std::endl;
            std::cout << "AFTER: this->unsorted_minimiser_it: " << (mixed_ptr{this->unsorted_minimiser_it} - mixed_ptr{window_values.data(), this}) << std::endl;
        }

        return minimiser_changed;
    }

protected:

    std::pair<value_type *, value_type *> _valid_sorted_data()
    {
        value_type * sorted_begin = this->window_values.data();
        return {sorted_begin, sorted_begin + window_size};
    }

    std::pair<value_type *, value_type *> _valid_unsorted_data()
    {
        value_type * unsorted_begin = this->window_values.data() + 1 + window_size;
        return {unsorted_begin, unsorted_begin + window_size};
    }

    void diagnostics()
    {
        std::cout << "stellar_minimiser_window:" << std::endl;
        std::cout << "\tsorted:   [";
        for (auto it = this->sorted_begin; it != this->sorted_end; ++it)
            std::cout << ((mixed_ptr)it == this->minimiser_it ? "*" : "") << *it << ", ";
        std::cout << "]" << std::endl;
        std::cout << "\tunsorted: [";
        for (auto it = this->unsorted_begin - 1; it != this->unsorted_end; ++it)
            std::cout << ((mixed_ptr)it == this->minimiser_it ? "*" : "") << (it == this->unsorted_minimiser_it ? "!" : "") << *it << ", ";
        std::cout << "]" << std::endl;
    }

    template <typename ptr_t>
    static constexpr ptr_t indexed_minimum(ptr_t const old_it, ptr_t const new_it)
    {
        return *new_it < *old_it ? new_it : old_it;
    }

    // Note recalculate_minimum will do:
    // In range (sorted_end, minimiser_it] recompute minimiser backwards (old minimiser would have left window)
    // In range (minimiser_it, sorted_begin] set old minimiser value as it is still active in that range
    void recalculate_minimum()
    {
        assert(mixed_ptr{this->sorted_end} + 1 == mixed_ptr{this->unsorted_begin});

        sorted_ptr sorted_it = this->sorted_end;
        unsorted_ptr unsorted_it = this->unsorted_end;
        unsorted_ptr unsorted_sentinel = this->unsorted_begin;

        // This is will be in the sorted memory region
        unsorted_ptr current_minimiser_it = this->unsorted_minimiser_it;

        // unsorted is non-empty
        assert(unsorted_it != unsorted_sentinel);

        // update pointer to be e in memory layout
        // S: sorted_it
        // U: unsorted_it
        // UE: unsorted_sentinel
        // UM: unsorted_minimiser_it
        // M: minimiser_it
        // [ s1, s2, ..sM.., sW, e*  , u1, u2, ..uM.., uW] (old memory)
        //                       S|UM, UE,       M       ] U (pointer)
        this->unsorted_minimiser_it = unsorted_sentinel - 1; // same as sorted_it
        this->minimiser_it = mixed_ptr{unsorted_sentinel - 1};
        // [ s1, s2, ..sM.., sW, uW  , u1, u2, ..uM.., uW] (new memory)
        *this->unsorted_minimiser_it = std::numeric_limits<value_type>::max();

        // construct minimiser from last to "first" element
        std::cout << "PRE_MINIMISER" << std::endl;
        for (; unsorted_it != current_minimiser_it; )
        {
            --unsorted_it;
            --sorted_it;
            std::cout << "U[" << (unsorted_it - unsorted_ptr{_valid_unsorted_data().first, this}) << "/" << (current_minimiser_it.ptr - _valid_unsorted_data().first) << "]: " << *unsorted_it << std::endl;
            std::cout << "S[" << (sorted_it - sorted_ptr{_valid_sorted_data().first, this}) << "] = " << *sorted_it << " := " << *unsorted_it << std::endl;
            *sorted_it = *unsorted_it;
            this->minimiser_it = indexed_minimum(this->minimiser_it, (mixed_ptr)sorted_it);
            *sorted_it = *this->minimiser_it;
            std::cout << "S[" << (sorted_it - sorted_ptr{_valid_sorted_data().first, this}) << "] = " << *sorted_it << " := " << *this->minimiser_it << std::endl;

            // vM+1, ..., vW: minimiser vi := min(ui, ..., uW)
            // [ s1, s2, ..., sM, sM+1..., sW, e*, u1, u2, ..., uM, uM+1 ..., uW] (old memory)
            //                    S            UM, UE,              U,          ] (pointer)
            // [ s1, s2, ..., sM, vM+1..., vW, uW, u1, u2, ..., uM, uM+1 ..., uW] (new memory)
        }

        std::cout << "POST_MINIMISER" << std::endl;
        // set old minimiser from its position until the "first" element
        for (; unsorted_it != unsorted_sentinel; )
        {
            --unsorted_it;
            --sorted_it;
            std::cout << "S[" << (sorted_it - sorted_ptr{_valid_sorted_data().first, this}) << "] = " << *sorted_it << " := " << *current_minimiser_it << std::endl;
            *sorted_it = *current_minimiser_it;
            mixed_ptr current_minimiser_sorted_it = mixed_ptr{current_minimiser_it} - window_size - 1;
            // assert(indexed_minimum((sorted_ptr)this->minimiser_it, sorted_it) == (sorted_ptr)current_minimiser_sorted_it);
            // this->minimiser_it = (mixed_ptr)current_minimiser_it;

            // vM+1, ..., vW: minimiser vi := min(ui, ..., uW)
            // [ s1, s2, ..., sM, sM+1..., sW, e*, u1,   u2, ..., uM, uM+1 ..., uW] (old memory)
            //   S                             UM, UE|U,                          ] (pointer)
            // [ uM, uM, ..., uM, vM+1..., vW, uW, u1,   u2, ..., uM, uM+1 ..., uW] (new memory)
        }

        // update pointer
        this->sorted_begin = sorted_it; // sorted region is now non-empty
        this->unsorted_end = this->unsorted_begin; // unsorted region is now empty
    }

    struct sorted_ptr;
    struct unsorted_ptr;

    template <typename derived_t>
    struct basic_ptr
    {
    private:
        friend derived_t;

        basic_ptr() = default;
        basic_ptr(basic_ptr &&) = default;
        basic_ptr(basic_ptr const &) = default;
        basic_ptr & operator=(basic_ptr &&) = default;
        basic_ptr & operator=(basic_ptr const &) = default;
    public:
        basic_ptr(value_type * ptr, stellar_minimiser_window * host_ptr)
            : ptr{ptr}, host_ptr{host_ptr}
        {
            _assert_bounds();
        }

        decltype(auto) operator*() { _assert_bounds(true); return *ptr; }
        decltype(auto) operator*() const { _assert_bounds(true); return *ptr; }

        void operator++() { ++ptr; _assert_bounds(); }
        void operator--() { --ptr; _assert_bounds(); }
        derived_t operator-(std::ptrdiff_t d) { return {ptr - d, host_ptr}; }
        derived_t operator+(std::ptrdiff_t d) { return {ptr + d, host_ptr}; }

        friend auto operator-(derived_t const & ptr1, derived_t const & ptr2)
        {
            return ptr1.ptr - ptr2.ptr;
        }

        friend bool operator==(derived_t const & ptr1, derived_t const & ptr2)
        {
            return ptr1.ptr == ptr2.ptr;
        }

        friend bool operator<(derived_t const & ptr1, derived_t const & ptr2)
        {
            return ptr1.ptr < ptr2.ptr;
        }

        void _assert_bounds(bool const dereferencable = false) const
        {
#ifndef NDEBUG
            auto in_range = [](value_type * ptr, value_type * begin, value_type * end)
            {
                return begin <= ptr && ptr <= end;
            };

            auto implies = [](bool const a, bool const b)
            {
                return !a || b;
            };

            value_type * window_begin = host_ptr->window_values.data();
            value_type * window_end = window_begin + host_ptr->window_values.size() + (dereferencable ? -1 : 0);
            {
                if (!in_range(ptr, window_begin, window_end))
                    std::cout << (dereferencable ? "*" : "") << "W: " << (window_begin - window_begin) << " <= " << (ptr - window_begin) << " < " << (window_end - window_begin) << std::endl;
                assert(in_range(ptr, window_begin, window_end));
            }

            {
                auto [sorted_begin, sorted_end] = host_ptr->_valid_sorted_data();
                sorted_end += (dereferencable ? -1 : 0);
                bool in_sorted = std::is_same_v<derived_t, sorted_ptr>;
                bool is_initialized =
                    in_range(sorted_begin, window_begin, window_end) &&
                    in_range(sorted_end, window_begin, window_end);
                if (!implies(in_sorted && is_initialized, in_range(ptr, sorted_begin, sorted_end)))
                    std::cout << (dereferencable ? "*" : "") << "S: " << (sorted_begin - window_begin) << " <= " << (ptr - window_begin) << " < " << (sorted_end - window_begin) << std::endl;
                assert(implies(in_sorted && is_initialized, in_range(ptr, sorted_begin, sorted_end)));
            }

            {
                auto [unsorted_begin, unsorted_end] = host_ptr->_valid_unsorted_data();
                --unsorted_begin;
                unsorted_end += (dereferencable ? -1 : 0);
                bool in_unsorted = std::is_same_v<derived_t, unsorted_ptr>;
                bool is_initialized =
                    in_range(unsorted_begin, window_begin, window_end) &&
                    in_range(unsorted_end, window_begin, window_end);
                if (!implies(in_unsorted && is_initialized, in_range(ptr, unsorted_begin, unsorted_end)))
                    std::cout << (dereferencable ? "*" : "") << "U: " << (unsorted_begin - window_begin) << " <= " << (ptr - window_begin) << " < " << (unsorted_end - window_begin) << std::endl;
                assert(implies(in_unsorted && is_initialized, in_range(ptr, unsorted_begin, unsorted_end)));
            }

            // std::cout << std::endl;
#endif // NDEBUG
        }

        value_type * ptr;

        stellar_minimiser_window * host_ptr;
    };

    struct sorted_ptr : public basic_ptr<sorted_ptr>
    {
        sorted_ptr() = default;
        sorted_ptr(value_type * ptr, stellar_minimiser_window * host_ptr)
            : basic_ptr<sorted_ptr>{ptr, host_ptr}
        {}
    };

    struct unsorted_ptr : public basic_ptr<unsorted_ptr>
    {
        unsorted_ptr() = default;
        unsorted_ptr(value_type * ptr, stellar_minimiser_window * host_ptr)
            : basic_ptr<unsorted_ptr>{ptr, host_ptr}
        {}
    };

    struct mixed_ptr : public basic_ptr<mixed_ptr>
    {
        mixed_ptr() = default;
        mixed_ptr(value_type * ptr, stellar_minimiser_window * host_ptr)
            : basic_ptr<mixed_ptr>{ptr, host_ptr}
        {}

        explicit mixed_ptr(sorted_ptr ptr)
            : basic_ptr<mixed_ptr>{ptr.ptr, ptr.host_ptr}
        {}
        explicit mixed_ptr(unsorted_ptr ptr)
            : basic_ptr<mixed_ptr>{ptr.ptr, ptr.host_ptr}
        {}

        explicit operator sorted_ptr() const
        {
            return sorted_ptr{this->ptr, this->host_ptr};
        }

        explicit operator unsorted_ptr() const
        {
            return unsorted_ptr{this->ptr, this->host_ptr};
        }
    };

    mixed_ptr minimiser_it{};
    unsorted_ptr unsorted_minimiser_it{};

    //!\brief Stored values per window. It is necessary to store them, because a shift can remove the current minimiser.
    std::vector<value_type> window_values{};

    sorted_ptr sorted_begin;
    sorted_ptr sorted_end;

    unsorted_ptr unsorted_begin;
    unsorted_ptr unsorted_end;

    size_t window_size;
};

} // namespace stellar::test
} // namespace stellar
