#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>
#include <stdexcept>

#ifndef NDEBUG
// #define STELLAR_MINIMISER_WINDOW_DEBUG
#endif

namespace stellar {
namespace test {

template <typename value_type>
struct stellar_minimiser_window
{
protected:
    struct sorted_ptr;
    struct unsorted_ptr;
    struct mixed_ptr;

public:
    // Memory Layout
    // s1, s2, ..., sW: sorted elements
    // e*: one element between sorted and unsorted
    // u1, u2, ..., uW: unsorted elements
    // [ s1, s2, ..., sW, e*, u1, u2, ..., uW]
    stellar_minimiser_window(size_t const window_size) :
        window_size{window_size}
    {
        window_values.resize(2 * window_size + 1);
        sorted_minimizer_stack.reserve(window_size);
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
        iterator_t new_input_sentinel = input_iterator + window_size;

        this->sorted_begin = {_valid_sorted_data().first, this};
        this->sorted_end = {std::copy(input_iterator, new_input_sentinel, this->sorted_begin.ptr), this};

        this->unsorted_begin = {_valid_unsorted_data().first, this};
        this->unsorted_end = this->unsorted_begin;

        *unsorted_infinity() = std::numeric_limits<value_type>::max();
        this->unsorted_minimiser_it = unsorted_infinity();

#ifdef STELLAR_MINIMISER_WINDOW_DEBUG
        std::cout << "INIT!" << std::endl;
#endif // STELLAR_MINIMISER_WINDOW_DEBUG
        recalculate_minimum(true);

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
#ifdef STELLAR_MINIMISER_WINDOW_DEBUG
        std::cout << "~~~~~~~~~~~~~~~~~~~" << std::endl;
        std::cout << "cyclic_push(" << new_value << "):" << std::endl;
#endif // STELLAR_MINIMISER_WINDOW_DEBUG
        ++this->sorted_begin;
        mixed_ptr previous_minimiser_it = this->minimiser_it;

        bool const sorted_range_empty = this->sorted_begin == this->sorted_end;
        bool const previous_minimiser_left_window = previous_minimiser_it < (mixed_ptr)this->sorted_begin;
        bool const previous_minimiser_was_in_sorted = previous_minimiser_it < (mixed_ptr)this->sorted_end;

        // append new element to unsorted list and update minimiser of the unsorted list
        *this->unsorted_end = new_value;
        if (!previous_minimiser_was_in_sorted)
        {
            // minimiser_in_unsorted: If the minimiser is already in the unsorted range, keep the current position
            // unless a REAL smaller value appears.
            // unsorted: [*5, 5, 6], 4, 4 <- should not change the position to the second 5
            //
            // ALSO
            //
            // unsorted: [5, 6, *4], 4 <- will just choose first 4 (as second one is not in window yet)
            // unsorted: [6, *4, 4] <- a consecutive 4 will not change minimiser position (as current number didn't left the window)
            this->unsorted_minimiser_it = indexed_minimum(this->unsorted_minimiser_it, this->unsorted_end);
        } else {
            // previous_minimiser_was_in_sorted: "consecutive" runs of the same value should be the last one [0, 1, 0, *0]
            this->unsorted_minimiser_it = indexed_minimum_less_equal(this->unsorted_minimiser_it, this->unsorted_end);
        }
        ++this->unsorted_end;

        // update overall minimiser
        if (sorted_range_empty)
        {
            this->minimiser_it = (mixed_ptr)this->unsorted_minimiser_it;
             // note: the minimiser can also not change, as unsorted_minimiser_it could have stayed the same, but will
             // be selected as overall minimiser, since no minimiser in sorted are left.
            bool minimiser_changed = previous_minimiser_it != this->minimiser_it;

#ifdef STELLAR_MINIMISER_WINDOW_DEBUG
            this->diagnostics();
            std::cout << "REBUILD!" << std::endl;
            std::cout << "BEFORE: this->unsorted_minimiser_it: " << mixed_ptr{this->unsorted_minimiser_it}._debug_position() << std::endl;
            std::cout << "BEFORE: this->minimiser_it: " << this->minimiser_it._debug_position() << std::endl;
            std::cout << "BEFORE: previous_minimiser_it: " << previous_minimiser_it._debug_position() << std::endl;
#endif // STELLAR_MINIMISER_WINDOW_DEBUG

            recalculate_minimum(false);

            return minimiser_changed;
        }
        else if (previous_minimiser_left_window)
        {
            assert(!sorted_minimizer_stack.empty());
            this->minimiser_it = (mixed_ptr)sorted_minimizer_stack.back();
            sorted_minimizer_stack.pop_back();

            // if sorted range and unsorted range have same "value", prefer the unsorted one
            // example:
            //  sorted_minimiser:   [[2]: 18, [1]: 8, ]
            //  sorted:   [*0, 8, 18, ]
            //  unsorted: [!2147483647, ]
            //
            //  cyclic_push(8)
            //  sorted_minimiser:   [[2]: 18, [1]: 8, ]
            //  sorted:   [*8, 18, ]
            //  unsorted: [2147483647, !8] <- prefer this eight, as 0 left the window and in that case we always
            //                                prefer the most right one.
            this->minimiser_it = indexed_minimum_less_equal(this->minimiser_it, (mixed_ptr)this->unsorted_minimiser_it);

#ifdef STELLAR_MINIMISER_WINDOW_DEBUG
            this->diagnostics();
            std::cout << "previous_minimiser_it: " << previous_minimiser_it._debug_position() << std::endl;
            std::cout << "this->minimiser_it: " << this->minimiser_it._debug_position() << std::endl;
#endif // STELLAR_MINIMISER_WINDOW_DEBUG

            return true;
        } else
        {
            // example:
            //  sorted:   [0, 0, 0, *0, ]
            //  unsorted: [!2147483647, ]
            //
            //  cyclic_push(0)
            //  sorted:   [0, 0, *0, ] <- the active minimiser should not be overwritten (since minimiser didn't leave the window)
            //  unsorted: [2147483647, !0] <- don't prefer this one
            this->minimiser_it = indexed_minimum(this->minimiser_it, (mixed_ptr)this->unsorted_minimiser_it);

#ifdef STELLAR_MINIMISER_WINDOW_DEBUG
            this->diagnostics();
            std::cout << "previous_minimiser_it: " << previous_minimiser_it._debug_position() << std::endl;
            std::cout << "this->minimiser_it: " << this->minimiser_it._debug_position() << std::endl;
#endif // STELLAR_MINIMISER_WINDOW_DEBUG

            bool minimiser_changed = previous_minimiser_it != this->minimiser_it;
            return minimiser_changed;
        }
    }

protected:

    unsorted_ptr unsorted_infinity()
    {
        return {_valid_unsorted_data().first - 1, this};
    }

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
#ifdef STELLAR_MINIMISER_WINDOW_DEBUG
        std::cout << "stellar_minimiser_window:" << std::endl;
        std::cout << "\tsorted_minimiser:   [";
        bool minimizer_seen = false;
        bool unsorted_minimizer_seen = false;
        for (sorted_ptr ptr : sorted_minimizer_stack)
        {
            std::cout << ((mixed_ptr)ptr == this->minimiser_it ? "*" : "") << "[" << ptr._debug_position() << "]: " << *ptr << ", ";
        }
        std::cout << "]" << std::endl;
        std::cout << "\tsorted:   [";
        for (auto it = this->sorted_begin; it != this->sorted_end; ++it)
        {
            minimizer_seen |= (mixed_ptr)it == this->minimiser_it;
            std::cout << ((mixed_ptr)it == this->minimiser_it ? "*" : "") << *it << ", ";
        }
        std::cout << "]" << std::endl;
        std::cout << "\tunsorted: [";
        for (auto it = this->unsorted_begin - 1; it != this->unsorted_end; ++it)
        {
            minimizer_seen |= (mixed_ptr)it == this->minimiser_it;
            unsorted_minimizer_seen |= it == this->unsorted_minimiser_it;
            std::cout << ((mixed_ptr)it == this->minimiser_it ? "*" : "") << (it == this->unsorted_minimiser_it ? "!" : "") << *it << ", ";
        }
        std::cout << "]" << std::endl;

        assert(minimizer_seen);
        assert(unsorted_minimizer_seen);
#endif // STELLAR_MINIMISER_WINDOW_DEBUG
    }

    template <typename ptr_t>
    static constexpr ptr_t indexed_minimum(ptr_t const old_it, ptr_t const new_it)
    {
        return *new_it < *old_it ? new_it : old_it;
    }

    template <typename ptr_t>
    static constexpr ptr_t indexed_minimum_less_equal(ptr_t const old_it, ptr_t const new_it)
    {
        return indexed_minimum(new_it, old_it);
    }

    // Note recalculate_minimum will do:
    // In range (sorted_end, minimiser_it] recompute minimiser backwards (old minimiser would have left window)
    // In range (minimiser_it, sorted_begin] set old minimiser value as it is still active in that range
    void recalculate_minimum(bool in_initialization)
    {
        sorted_minimizer_stack.clear();

        assert(mixed_ptr{this->sorted_end} + 1 == mixed_ptr{this->unsorted_begin});

        // unsorted_ptr const current_minimiser_unsorted_it = this->unsorted_minimiser_it;
        std::ptrdiff_t minimiser_position = this->unsorted_minimiser_it - this->unsorted_begin;
#ifdef STELLAR_MINIMISER_WINDOW_DEBUG
        std::cout << "PRE_MINIMISER" << std::endl;
#endif // STELLAR_MINIMISER_WINDOW_DEBUG

        // copy unsorted elements into sorted elements
        if (!in_initialization)
        {
            // unsorted is non-empty
            assert(this->unsorted_begin != this->unsorted_end);

            // there exists a unsorted minimiser
            assert(!(this->unsorted_minimiser_it < this->unsorted_begin));

            sorted_ptr sorted_begin{_valid_sorted_data().first, this};
            sorted_ptr current_minimiser_sorted_it = sorted_begin + minimiser_position;
            // there are exactly "window_size" many elements in the unsorted list.
            // Note: Actually there can be less elements in sorted, if the initial range is to small, but this can only
            // happen if in_initialization == true.
            //
            // in_initialization == false means that we trigger a rebuild by cyclic_push and that means that enough
            // initial elements where in the queue.
            assert(static_cast<size_t>(this->unsorted_end - this->unsorted_begin) == window_size);
            assert(sorted_begin == this->sorted_end - (this->unsorted_end - this->unsorted_begin));
            sorted_ptr sorted_end
            {
                std::copy(this->unsorted_minimiser_it.ptr, this->unsorted_end.ptr, current_minimiser_sorted_it.ptr),
                this
            };
            assert(this->sorted_end == sorted_end);

            // update pointer
            this->sorted_begin = sorted_begin; // sorted region is now non-empty
            this->unsorted_end = this->unsorted_begin; // unsorted region is now empty
            this->unsorted_minimiser_it = unsorted_infinity();
        }

        // construct minimiser from last to "first" element
        {
            mixed_ptr minimiser_it = mixed_ptr{unsorted_infinity()};
            sorted_ptr current_minimiser_sorted_it = this->sorted_begin + (in_initialization ? 0 : minimiser_position);
            assert(current_minimiser_sorted_it != this->sorted_end); // range is not empty
            for (sorted_ptr sorted_it = this->sorted_end; sorted_it != current_minimiser_sorted_it; )
            {
                --sorted_it;
                sorted_ptr new_minimiser_it = (sorted_ptr)indexed_minimum(minimiser_it, (mixed_ptr)sorted_it);

#ifdef STELLAR_MINIMISER_WINDOW_DEBUG
                std::cout << "S[" << sorted_it._debug_position() << "] = " << *sorted_it << std::endl;
                std::cout << "new_minimiser_it: S[" << new_minimiser_it._debug_position() << "]: " << *new_minimiser_it << std::endl;
                if ((mixed_ptr)new_minimiser_it != minimiser_it)
                    std::cout << "minimiser changed" << std::endl;
#endif // STELLAR_MINIMISER_WINDOW_DEBUG

                if ((mixed_ptr)new_minimiser_it != minimiser_it)
                    sorted_minimizer_stack.push_back(new_minimiser_it);

                minimiser_it = (mixed_ptr)new_minimiser_it;
            }

            // The old minimiser must be re-added if it is different from recomputed one.
            //
            // example:
            // before:
            //    unsorted [5, !3, 4, 3, 3, 3]
            //
            // after:
            //    sorted: [5, 3, 4, 3, 3, *3] <- overall minimizer is the rightest one
            //    sorted: [5, *3, 4, 3, 3, *3] <- but as the old minimiser didn't leave the window yet, we can't ignore it
            if (!in_initialization && (mixed_ptr)current_minimiser_sorted_it != minimiser_it)
                sorted_minimizer_stack.push_back(current_minimiser_sorted_it);

            this->minimiser_it = (mixed_ptr)sorted_minimizer_stack.back();
            sorted_minimizer_stack.pop_back();
        }

#ifdef STELLAR_MINIMISER_WINDOW_DEBUG
        std::cout << "this->minimiser_it: S[" << this->minimiser_it._debug_position() << "]: " << *this->minimiser_it << std::endl;
        std::cout << "POST_MINIMISER" << std::endl;
#endif // STELLAR_MINIMISER_WINDOW_DEBUG
    }

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
        basic_ptr(value_type * ptr, [[maybe_unused]] stellar_minimiser_window const * host_ptr) :
            ptr{ptr}
#ifndef NDEBUG
            , host_ptr{host_ptr}
#endif // NDEBUG
        {
            _assert_bounds();
        }

        decltype(auto) operator*() { _assert_bounds(true); return *ptr; }
        decltype(auto) operator*() const { _assert_bounds(true); return *ptr; }

        void operator++() { ++ptr; _assert_bounds(); }
        void operator--() { --ptr; _assert_bounds(); }
        derived_t operator-(std::ptrdiff_t d) { return {ptr - d, _host_ptr()}; }
        derived_t operator+(std::ptrdiff_t d) { return {ptr + d, _host_ptr()}; }

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

        stellar_minimiser_window const * _host_ptr() const
        {
#ifndef NDEBUG
            return host_ptr;
#endif // NDEBUG
            return nullptr;
        }

        std::ptrdiff_t _debug_position() const
        {
#ifndef NDEBUG
            auto [sorted_begin, sorted_end] = host_ptr->_valid_sorted_data();
            if (std::is_same_v<derived_t, sorted_ptr> && sorted_begin <= ptr && ptr <= sorted_end)
                return ptr - sorted_begin;

            auto [unsorted_begin, unsorted_end] = host_ptr->_valid_unsorted_data();
            if (std::is_same_v<derived_t, unsorted_ptr> && unsorted_begin - 1 <= ptr && ptr <= unsorted_end)
                return ptr - unsorted_begin;

            value_type * window_begin = host_ptr->window_values.data();
            value_type * window_end = window_begin + host_ptr->window_values.size();
            if (window_begin <= ptr && ptr <= window_end)
                return ptr - window_begin;

#endif
            throw std::out_of_range{"_debug_position called, but ptr is invalid"};
        }

        void _assert_bounds([[maybe_unused]] bool const dereferencable = false) const
        {
#ifdef STELLAR_MINIMISER_WINDOW_DEBUG
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
                assert(is_initialized);
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
                assert(is_initialized);
                if (!implies(in_unsorted && is_initialized, in_range(ptr, unsorted_begin, unsorted_end)))
                    std::cout << (dereferencable ? "*" : "") << "U: " << (unsorted_begin - window_begin) << " <= " << (ptr - window_begin) << " < " << (unsorted_end - window_begin) << std::endl;
                assert(implies(in_unsorted && is_initialized, in_range(ptr, unsorted_begin, unsorted_end)));
            }

            // std::cout << std::endl;
#endif // STELLAR_MINIMISER_WINDOW_DEBUG
        }

        value_type * ptr;

#ifndef NDEBUG
        stellar_minimiser_window const * host_ptr{nullptr};
#endif // NDEBUG
    };

    struct sorted_ptr : public basic_ptr<sorted_ptr>
    {
        sorted_ptr() = default;
        sorted_ptr(value_type * ptr, stellar_minimiser_window const * host_ptr)
            : basic_ptr<sorted_ptr>{ptr, host_ptr}
        {}
    };

    struct unsorted_ptr : public basic_ptr<unsorted_ptr>
    {
        unsorted_ptr() = default;
        unsorted_ptr(value_type * ptr, stellar_minimiser_window const * host_ptr)
            : basic_ptr<unsorted_ptr>{ptr, host_ptr}
        {}
    };

    struct mixed_ptr : public basic_ptr<mixed_ptr>
    {
        mixed_ptr() = default;
        mixed_ptr(value_type * ptr, stellar_minimiser_window const * host_ptr)
            : basic_ptr<mixed_ptr>{ptr, host_ptr}
        {}

        explicit mixed_ptr(sorted_ptr ptr)
            : basic_ptr<mixed_ptr>{ptr.ptr, ptr._host_ptr()}
        {}
        explicit mixed_ptr(unsorted_ptr ptr)
            : basic_ptr<mixed_ptr>{ptr.ptr, ptr._host_ptr()}
        {}

        explicit operator sorted_ptr() const
        {
            return sorted_ptr{this->ptr, this->_host_ptr()};
        }

        explicit operator unsorted_ptr() const
        {
            return unsorted_ptr{this->ptr, this->_host_ptr()};
        }
    };

    mixed_ptr minimiser_it{};
    unsorted_ptr unsorted_minimiser_it{};

    std::vector<sorted_ptr> sorted_minimizer_stack{};

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
