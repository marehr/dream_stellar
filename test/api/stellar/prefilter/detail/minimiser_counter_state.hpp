
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
