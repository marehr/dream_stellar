
#pragma once

#include <stellar/prefilter/detail/valik_local_query_prefilter_index.hpp>

template <typename TAlphabet>
struct minimiser_counter_state
{
    using prefilter_index_t = stellar::detail::valik_local_query_prefilter_index;
    using ibf_t = decltype(prefilter_index_t{}._ibf);
    using counting_agent_t = typename ibf_t::counting_agent_type<uint16_t>;
    using counting_vector_t = seqan3::counting_vector<uint16_t>;

    void init() {};

    bool next()
    {
        if ((_current_position + subquery_size) > length(query_sequence))
            return false;

        size_t new_current_position = _current_position;
        std::vector<size_t> new_minimiser_hashes{};

        while ((_current_position + subquery_size) <= length(query_sequence))
        {
            auto subquery_sequence = subquery();
            auto minimiser_hashes_view = lqpindex.hash_view(subquery_sequence) | std::views::common;
            new_minimiser_hashes = std::vector<size_t>(minimiser_hashes_view.begin(), minimiser_hashes_view.end());

            if (!std::ranges::equal(_minimiser_hashes, new_minimiser_hashes))
                break;

            ++_current_position;
        }

        _next_unique_position = _current_position;
        _current_position = new_current_position;

        if (new_current_position == _next_unique_position)
            _minimiser_hashes = std::move(new_minimiser_hashes);

        [[maybe_unused]] auto const & ref = counting_agent.bulk_count(_minimiser_hashes);

        ++_current_position;

        return true;
    }

    size_t next_unique_position() const
    {
        return _next_unique_position;
    }

    auto subquery() const
    {
        return seqan::infix(query_sequence, _current_position, _current_position + subquery_size);
    }

    std::vector<size_t> const & minimiser_hashes() const &
    {
        return _minimiser_hashes;
    }

    size_t minimiser_count() const
    {
        return _minimiser_hashes.size();
    }

    counting_vector_t const & bin_counts() const &
    {
        return counting_agent.result_buffer;
    }

    stellar::detail::valik_local_query_prefilter_index const & lqpindex;
    seqan::String<TAlphabet> const & query_sequence;
    std::size_t subquery_size = 0;

    std::size_t _next_unique_position = 0;
    std::size_t _current_position = 0;
    std::vector<size_t> _minimiser_hashes{};
    counting_agent_t counting_agent = lqpindex._ibf.counting_agent();
};

template <typename TAlphabet>
struct minimiser_counter_state2
{
    using prefilter_index_t = stellar::detail::valik_local_query_prefilter_index;
    using ibf_t = decltype(prefilter_index_t{}._ibf);
    using counting_agent_t = typename ibf_t::counting_agent_type<uint16_t>;
    using counting_vector_t = seqan3::counting_vector<uint16_t>;

    void init()
    {
        auto minimiser_hashes_view = lqpindex.hash_view(query_sequence);
        auto it = minimiser_hashes_view.begin();
        auto it_base = it.base();
        auto sentinel = minimiser_hashes_view.end();

        size_t window_size = lqpindex._window_size.get();
        size_t window_subquery_size = subquery_size - window_size;

        for (; it != sentinel; ++it)
        {
            // subquery_size: 15 (= min_length)
            // kmer_size: 4
            // kmer_subquery_size: 12 = 15 - 4 + 1
            // window_subquery_size: 9 = 15 - 6
            // chars:  [0|1|2|3|4|5|6|7|8|9|10|11|12|13|14]
            // kmer:   [   0   |1|2|3|4|5|6| 7| 8| 9|10|11]
            // window: [     0     |1|2|3|4| 5| 6| 7| 8| 9]
            auto kmer_begin_position = it.base() - it_base;
            // auto end_position = kmer_begin_position + window_size;
            _complete_minimiser_kmer_begin_positions.push_back(kmer_begin_position);
            // _complete_minimiser_end_positions.push_back(end_position);
            _complete_minimiser_hashes.push_back(*it);
        }

        seqan3::debug_stream << "init::_complete_minimiser_hashes: " << _complete_minimiser_hashes << std::endl;
        seqan3::debug_stream << "init::_complete_minimiser_kmer_begin_positions: " << _complete_minimiser_kmer_begin_positions << std::endl;

        auto upper_bound_it = std::upper_bound(
            _complete_minimiser_kmer_begin_positions.begin(),
            _complete_minimiser_kmer_begin_positions.end(),
            window_subquery_size);

        size_t upper_bound_distance = upper_bound_it - _complete_minimiser_kmer_begin_positions.begin();

        // _update_minimiser_begin_end_position
        _minimiser_begin_position = 0u;
        _minimiser_end_position = upper_bound_distance;

        seqan3::debug_stream << "init::subquery_size: " << subquery_size << std::endl;
        seqan3::debug_stream << "init::window_subquery_size: " << window_subquery_size << std::endl;
        seqan3::debug_stream << "init::_minimiser_begin_position: " << _minimiser_begin_position << std::endl;
        seqan3::debug_stream << "init::_minimiser_end_position: " << _minimiser_end_position << std::endl;

        size_t upper_bound_value =
            (upper_bound_it == _complete_minimiser_kmer_begin_positions.end()) ?
            std::numeric_limits<size_t>::max() :
            *upper_bound_it;
        std::cout << "upper_bound(" << window_subquery_size << "): < value: " << (upper_bound_value) << ", distance: " << upper_bound_distance << std::endl;

        _update_counting_vector();
        _update_next_event_in();
    }

    void _update_counting_vector()
    {
        [[maybe_unused]] auto const & ref = counting_agent.bulk_count(_minimiser_hashes_span());
    }

    bool at_end()
    {
        return (_current_position + subquery_size) > length(query_sequence);
    }

    void _update_minimiser_begin_end_position()
    {
        size_t next_event_in = std::min(_next_begin_event_in, _next_end_event_in);
        bool is_begin_event = (next_event_in == _next_begin_event_in);
        bool is_end_event = (next_event_in == _next_end_event_in);

        seqan3::debug_stream << "is_begin_event: " << (is_begin_event ? "true" : "false") << std::endl;
        seqan3::debug_stream << "is_end_event: " << (is_end_event ? "true" : "false") << std::endl;

        seqan3::debug_stream << "_minimiser_begin_position::before: " << _minimiser_begin_position << std::endl;
        seqan3::debug_stream << "_minimiser_end_position::before: " << _minimiser_end_position << std::endl;

        if (is_begin_event)
        {
            ++_minimiser_begin_position;
        }

        if (is_end_event)
        {
            ++_minimiser_end_position;
        }
        seqan3::debug_stream << "_minimiser_begin_position: " << _minimiser_begin_position << std::endl;
        seqan3::debug_stream << "_minimiser_end_position: " << _minimiser_end_position << std::endl;
    }

    void _update_next_event_in()
    {
        size_t window_size = lqpindex._window_size.get();

        seqan3::debug_stream << "_update_next_event_in::_current_position: " << _current_position << std::endl;
        size_t current_begin_position = _current_position;
        size_t current_end_position = _current_position + subquery_size;
        size_t next_begin_position = _complete_minimiser_kmer_begin_positions[_minimiser_begin_position + 1];
        size_t next_end_position = _complete_minimiser_kmer_begin_positions[_minimiser_end_position] + window_size;

        _next_begin_event_in = next_begin_position - current_begin_position;
        _next_end_event_in = next_end_position - current_end_position;

        size_t next_event_in = std::min(_next_begin_event_in, _next_end_event_in);

        seqan3::debug_stream << "_update_next_event_in::current_begin_position: " << current_begin_position << std::endl;
        seqan3::debug_stream << "_update_next_event_in::current_end_position: " << current_end_position << std::endl;
        seqan3::debug_stream << "_update_next_event_in::next_begin_position: " << next_begin_position << std::endl;
        seqan3::debug_stream << "_update_next_event_in::next_end_position: " << next_end_position << std::endl;
        seqan3::debug_stream << "_update_next_event_in::_next_begin_event_in: " << _next_begin_event_in << std::endl;
        seqan3::debug_stream << "_update_next_event_in::_next_end_event_in: " << _next_end_event_in << std::endl;
        seqan3::debug_stream << "_update_next_event_in::next_event_in: " << next_event_in << std::endl;

        seqan3::debug_stream << "_update_next_event_in::_next_unique_position::before: " << _next_unique_position << std::endl;
        _next_unique_position += next_event_in;
        seqan3::debug_stream << "_update_next_event_in::_next_unique_position: " << _next_unique_position << std::endl;
    }

    void next()
    {
        if (at_end())
            return;

        ++_current_position;

        seqan3::debug_stream << std::endl;
        seqan3::debug_stream << "_current_position: " << _current_position << std::endl;
        seqan3::debug_stream << "_next_unique_position: " << _next_unique_position << std::endl;

        bool update_minimiser_hashes = _next_unique_position == _current_position;
        seqan3::debug_stream << "update_minimiser_hashes: " << (update_minimiser_hashes ? "true" : "false") << std::endl;

        if (update_minimiser_hashes)
        {
            seqan3::debug_stream << "_complete_minimiser_hashes: " << _complete_minimiser_hashes << std::endl;
            seqan3::debug_stream << "_complete_minimiser_kmer_begin_positions: " << _complete_minimiser_kmer_begin_positions << std::endl;

            _update_minimiser_begin_end_position();
            _update_counting_vector();
            _update_next_event_in();
        }
    }

    size_t next_unique_position() const
    {
        return _next_unique_position;
    }

    size_t position() const
    {
        return _current_position;
    }

    auto subquery() const
    {
        return seqan::infix(query_sequence, _current_position, _current_position + subquery_size);
    }

    std::span<size_t const> _minimiser_hashes_span() const
    {
        return std::span<size_t const>{_complete_minimiser_hashes}
            .subspan(_minimiser_begin_position, minimiser_count());
    }

    std::vector<size_t> minimiser_hashes()
    {
        std::span<size_t const> minimiser_hashes_span = _minimiser_hashes_span();
        return std::vector<size_t>(minimiser_hashes_span.begin(), minimiser_hashes_span.end());
    }

    size_t minimiser_count() const
    {
        return _minimiser_end_position - _minimiser_begin_position;
    }

    counting_vector_t const & bin_counts() const &
    {
        return counting_agent.result_buffer;
    }

    stellar::detail::valik_local_query_prefilter_index const & lqpindex;
    seqan::String<TAlphabet> const & query_sequence;
    std::size_t subquery_size = 0;

    std::size_t _current_position = 0;
    std::size_t _next_unique_position = 0;

    std::size_t _minimiser_begin_position{};
    std::size_t _minimiser_end_position{};

    std::size_t _next_begin_event_in = 0;
    std::size_t _next_end_event_in = 0;

    // std::vector<size_t> _minimiser_hashes{}; // not needed

    std::vector<size_t> _complete_minimiser_hashes{};
    std::vector<size_t> _complete_minimiser_kmer_begin_positions{};
    counting_agent_t counting_agent = lqpindex._ibf.counting_agent();
};
