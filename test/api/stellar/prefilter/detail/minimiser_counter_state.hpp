
#pragma once

#include <stellar/prefilter/detail/valik_local_query_prefilter_index.hpp>

template <typename TAlphabet>
struct minimiser_counter_state
{
    using prefilter_index_t = stellar::detail::valik_local_query_prefilter_index;
    using ibf_t = decltype(prefilter_index_t{}._ibf);
    using counting_agent_t = typename ibf_t::counting_agent_type<uint16_t>;
    using counting_vector_t = seqan3::counting_vector<uint16_t>;

    bool next()
    {
        if ((_current_position + subquery_size) > length(query_sequence))
            return false;

        auto subquery_sequence = subquery();

        auto minimiser_hashes_view = lqpindex.hash_view(subquery_sequence) | std::views::common;
        _minimiser_hashes = std::vector<size_t>(minimiser_hashes_view.begin(), minimiser_hashes_view.end());
        _minimiser_count = _minimiser_hashes.size();

        [[maybe_unused]] auto const & ref = counting_agent.bulk_count(_minimiser_hashes);

        ++_current_position;

        return true;
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
        return _minimiser_count;
    }

    counting_vector_t const & bin_counts() const &
    {
        return counting_agent.result_buffer;
    }

    stellar::detail::valik_local_query_prefilter_index const & lqpindex;
    seqan::String<TAlphabet> const & query_sequence;
    std::size_t subquery_size = 0;

    std::size_t _current_position = 0;
    std::size_t _minimiser_count = 0;
    std::vector<size_t> _minimiser_hashes{};
    counting_agent_t counting_agent = lqpindex._ibf.counting_agent();
};
