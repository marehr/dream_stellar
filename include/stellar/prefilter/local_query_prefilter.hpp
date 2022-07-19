
#pragma once

#include <seqan3/std/span>

#include <stellar/stellar_database_segment.hpp>

#include <stellar/prefilter/detail/valik_local_query_prefilter_index.hpp>
#include <stellar/prefilter/detail/valik_local_prefilter.hpp>
// #include <stellar/prefilter/database_agent_distributer.hpp>
#include <stellar/prefilter/prefilter.hpp>
#include <stellar/prefilter/prefilter_agent.hpp>
#include <stellar/prefilter/prefilter_agents.hpp>

namespace stellar
{

// Prefilter: filters out the queries for each database sequence.
// In this case it does not filter out any queries and only builds a single query filter
template <typename TAlphabet>
struct LocalQueryPrefilter : stellar::prefilter<TAlphabet>
{
    using TQueryFilter = stellar::StellarSwiftPattern<TAlphabet>;

    LocalQueryPrefilter() = delete;

    LocalQueryPrefilter(
        stellar::IndexOptions const options,
        StringSet<String<TAlphabet> > const & databases,
        StringSet<String<TAlphabet> > const & queries,
        stellar::detail::valik_local_query_prefilter_index queryPrefilterIndex
    ) :
        _options{options},
        _databases{databases},
        _queries{queries},
        _queryPrefilterIndex{std::move(queryPrefilterIndex)}
    {}

    using TQuerySegment = seqan::Segment<seqan::String<TAlphabet> const, seqan::InfixSegment>;
    using TQuerySegments = std::span<TQuerySegment const>;

    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    using TDatabaseSegments = std::span<TDatabaseSegment const>;

    class Agent;

    virtual stellar::prefilter_agents<TAlphabet> agents(size_t const agentCount, stellar::EPSMatchOptions const epsMatchOptions) override
    {
        assert(_querySegments.size() == 0u); // calling agents twice is forbidden

        // _databaseSegments = _splitter.split(_databases, minLength);
        // TODO: use SPLITTER
        for (String<TAlphabet> const & query: _queries)
            if (length(query) >= epsMatchOptions.minLength)
                _querySegments.emplace_back(query, 0u, length(query));

        stellar::prefilter_agents<TAlphabet> agents{};
        agents.reserve(agentCount);

        // TODO: use DISTRIBUTER
        TQuerySegments allQuerySegments{_querySegments.begin(), _querySegments.end()};
        agents.template emplace_back<Agent>(*this, std::move(allQuerySegments), epsMatchOptions.minLength);
        // for (TDatabaseSegments & segments: DatabaseAgentDistributer::distribute(_databaseSegments, agentCount))
        //     agents.emplace_back(*this, std::move(segments));

        return agents;
    }

private:
    stellar::IndexOptions const _options{};
    StringSet<String<TAlphabet> > const & _databases{};
    StringSet<String<TAlphabet> > const & _queries{};
    std::vector<TQuerySegment> _querySegments{};
    stellar::detail::valik_local_query_prefilter_index _queryPrefilterIndex{};
};

template <typename TAlphabet>
struct LocalQueryPrefilter<TAlphabet>::Agent : prefilter_agent<TAlphabet>
{
    friend LocalQueryPrefilter;

    Agent(LocalQueryPrefilter & prefilter, TQuerySegments querySegments, size_t const minLength)
        : _prefilter_ptr{&prefilter}, _querySegments{std::move(querySegments)}, _minLength{minLength}
    {}

    using TBase = prefilter_agent<TAlphabet>;
public:

    using TPrefilterCallback = typename TBase::TPrefilterCallback;

    virtual void prefilter(TPrefilterCallback callback) override
    {
        // TODO: we somehow need to calculate this
        // TODO: not sure if this expression epsilon * minLength is correct
        size_t const maxError = 1u;
        assert(_prefilter_ptr != nullptr);

        auto binIDToDatabaseSegment = [&](size_t const binID) -> TDatabaseSegment
        {
            size_t const databaseID = binID; // TODO, retrieve translation from index structure
            auto const & database = _prefilter_ptr->_databases[databaseID];
            return TDatabaseSegment{database, 0u, length(database)};
        };

        for (TQuerySegment querySegment : _querySegments)
        {
            seqan::String<TAlphabet> const & query = host(querySegment);
            std::vector<TDatabaseSegment> _databaseSegments{};

            // This was used to get a std::cout if a bin_id met a threshold
#define STELLAR_DEBUG_PREFILTER_ALWAYS_MATCH 0

#if STELLAR_DEBUG_PREFILTER_ALWAYS_MATCH
            for (seqan::String<TAlphabet> const & database: _prefilter_ptr->_databases)
                _databaseSegments.emplace_back(database, 0u, length(database));
#endif // STELLAR_DEBUG_PREFILTER_ALWAYS_MATCH

#if STELLAR_DEBUG_LOCAL_PREFILTER
            std::cout << "valik_local_prefilter" << std::endl;
#endif // STELLAR_DEBUG_LOCAL_PREFILTER
            stellar::detail::valik_local_prefilter(
                query, // TODO: pass querySegment directly
                _prefilter_ptr->_queryPrefilterIndex,
                _minLength,
                maxError,
                [&](size_t const bin_id) {
#if STELLAR_DEBUG_LOCAL_PREFILTER
                    std::cout << "bin_id(" << bin_id << ") threshold met" << std::endl;
#endif // STELLAR_DEBUG_LOCAL_PREFILTER
#if !STELLAR_DEBUG_PREFILTER_ALWAYS_MATCH
                    _databaseSegments.emplace_back(binIDToDatabaseSegment(bin_id));
#endif // !STELLAR_DEBUG_PREFILTER_ALWAYS_MATCH
                });

            if (_databaseSegments.size() > 0)
            {
                std::span<TQuerySegment> querySegments{&querySegment, 1};
                stellar::StellarIndex<TAlphabet> index{querySegments, _prefilter_ptr->_options};
                // TQueryFilter _queryFilter{querySegment};
                auto pattern = index.createSwiftPattern();
                callback(_databaseSegments, pattern);
            }
        }
    }

private:
    LocalQueryPrefilter * _prefilter_ptr{nullptr};
    TQuerySegments _querySegments{};
    size_t _minLength{};
};

} // namespace stellar
