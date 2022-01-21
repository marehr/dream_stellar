
#pragma once

#include <seqan3/std/span>

#include <stellar/stellar_types.hpp>

#include <stellar/prefilter/database_agent_distributer.hpp>

namespace stellar
{

// Prefilter: filters out the queries for each database sequence.
// In this case it does not filter out any queries and only builds a single query filter
template <typename TAlphabet, typename TQueryFilter_, typename TAgentSplitter>
struct NoQueryPrefilter
{
    using TQueryFilter = TQueryFilter_;

    NoQueryPrefilter() = delete;

    NoQueryPrefilter(
        StringSet<String<TAlphabet> > const & databases,
        TQueryFilter queryFilter,
        TAgentSplitter agentSplitter = TAgentSplitter{}
    ) :
        _databases{databases},
        _queryFilter{std::move(queryFilter)},
        _splitter{std::move(agentSplitter)}
    {}

    template <typename TStringSet = StringSet<String<TAlphabet> > >
    NoQueryPrefilter(
        StringSet<String<TAlphabet> > const & databases,
        std::enable_if_t<!std::is_same_v<TQueryFilter, TStringSet const &>, TStringSet const &> queries
    ) :
        NoQueryPrefilter{databases, TQueryFilter{queries}}
    {}

    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    using TDatabaseSegments = std::span<TDatabaseSegment const>;

    class Agent;

    std::vector<Agent> agents(size_t const agentCount, size_t const minLength)
    {
        assert(_databaseSegments.size() == 0u); // calling agents twice is forbidden

        _databaseSegments = _splitter.split(_databases, minLength);

        std::vector<Agent> agents{};
        agents.reserve(agentCount);

        for (TDatabaseSegments & segments: DatabaseAgentDistributer::distribute(_databaseSegments, agentCount))
            agents.emplace_back(*this, std::move(segments));

        return agents;
    }

private:
    StringSet<String<TAlphabet> > const & _databases{};
    std::vector<TDatabaseSegment> _databaseSegments{};
    TQueryFilter _queryFilter{};
    TAgentSplitter _splitter{};
};

template <typename TAlphabet, typename TQueryFilter, typename TAgentSplitter>
struct NoQueryPrefilter<TAlphabet, TQueryFilter, TAgentSplitter>::Agent
{
    friend NoQueryPrefilter;

    Agent(NoQueryPrefilter & prefilter, TDatabaseSegments databaseSegments)
        : _prefilter_ptr{&prefilter}, _databaseSegments{std::move(databaseSegments)}
    {}
public:

    template <typename functor_t>
    void prefilter(functor_t && functor)
    {
        functor(_databaseSegments, _prefilter_ptr->_queryFilter);
    }

private:
    NoQueryPrefilter * _prefilter_ptr;
    TDatabaseSegments _databaseSegments{};
};


} // namespace stellar
