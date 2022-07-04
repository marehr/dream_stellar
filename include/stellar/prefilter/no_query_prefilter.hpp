
#pragma once

#include <seqan3/std/span>
#include <variant>

#include <stellar/options/index_options.hpp>

#include <stellar/prefilter/database_agent_distributer.hpp>
#include <stellar/prefilter/prefilter.hpp>
#include <stellar/prefilter/prefilter_agent.hpp>
#include <stellar/prefilter/prefilter_agents.hpp>

namespace stellar
{

// Prefilter: filters out the queries for each database sequence.
// In this case it does not filter out any queries and only builds a single query filter
template <typename TAlphabet, typename TAgentSplitter>
struct NoQueryPrefilter : stellar::prefilter<TAlphabet>
{
    using TQueryFilter = stellar::StellarSwiftPattern<TAlphabet>;

    NoQueryPrefilter() = delete;

    NoQueryPrefilter(
        StringSet<String<TAlphabet> > const & databases,
        TQueryFilter queryFilter,
        TAgentSplitter agentSplitter = TAgentSplitter{}
    ) :
        _databases{databases},
        _queryQGramIndex{},
        _queryFilter{std::move(queryFilter)}, // stellar::StellarSwiftPattern
        _splitter{std::move(agentSplitter)}
    {}

    template <typename TStringSet = StringSet<String<TAlphabet> > >
    NoQueryPrefilter(
        stellar::IndexOptions const options,
        StringSet<String<TAlphabet> > const & databases,
        std::enable_if_t<!std::is_same_v<TQueryFilter, TStringSet const &>, TStringSet const &> queries,
        TAgentSplitter agentSplitter = TAgentSplitter{}
    ) :
        _databases{databases},
        _queryQGramIndex{queries, options},
        _queryFilter{_queryQGramIndex.createSwiftPattern()},
        _splitter{std::move(agentSplitter)}
    {}

    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    using TDatabaseSegments = std::span<TDatabaseSegment const>;

    class Agent;

    virtual stellar::prefilter_agents<TAlphabet> agents(size_t const agentCount, stellar::EPSMatchOptions const epsMatchOptions) override
    {
        assert(_databaseSegments.size() == 0u); // calling agents twice is forbidden

        _databaseSegments = _splitter.split(_databases, epsMatchOptions.minLength);

        stellar::prefilter_agents<TAlphabet> agents{};
        agents.reserve(agentCount);

        for (TDatabaseSegments & segments: DatabaseAgentDistributer::distribute(_databaseSegments, agentCount))
            agents.template emplace_back<Agent>(*this, std::move(segments));

        return agents;
    }

private:

    struct QGramIndexVariant : std::variant<stellar::StellarQGramIndex<TAlphabet> *, stellar::StellarIndex<TAlphabet>>
    {
        using TBase = std::variant<stellar::StellarQGramIndex<TAlphabet> *, stellar::StellarIndex<TAlphabet>>;

        QGramIndexVariant() = default;

        QGramIndexVariant(StringSet<String<TAlphabet> > const & queries, stellar::IndexOptions const & options)
            : TBase{}
        {
            this->template emplace<stellar::StellarIndex<TAlphabet>>(queries, options);
        }

        QGramIndexVariant(stellar::StellarQGramIndex<TAlphabet> & qgramIndex)
            : TBase{&qgramIndex}
        {}

        stellar::StellarSwiftPattern<TAlphabet> createSwiftPattern()
        {
            return std::visit([](auto && index) -> stellar::StellarSwiftPattern<TAlphabet>
            {
                if constexpr (std::is_same_v<decltype(index), stellar::StellarQGramIndex<TAlphabet> * &>)
                    return stellar::StellarSwiftPattern<TAlphabet>{*index};
                else //if constexpr (std::is_same_v<decltype(index), stellar::StellarIndex<TAlphabet> &>)
                    return index.createSwiftPattern();
            }, *static_cast<TBase *>(this));
        }
    };

    StringSet<String<TAlphabet> > const & _databases;
    QGramIndexVariant _queryQGramIndex;
    std::vector<TDatabaseSegment> _databaseSegments{};
    TQueryFilter _queryFilter{};
    TAgentSplitter _splitter{};
};

template <typename TAlphabet, typename TAgentSplitter>
struct NoQueryPrefilter<TAlphabet, TAgentSplitter>::Agent : public prefilter_agent<TAlphabet>
{
    friend NoQueryPrefilter;

    Agent(NoQueryPrefilter & prefilter, TDatabaseSegments databaseSegments)
        : _prefilter_ptr{&prefilter}, _databaseSegments{std::move(databaseSegments)}
    {}

    using TBase = prefilter_agent<TAlphabet>;
public:

    using TPrefilterCallback = typename TBase::TPrefilterCallback;

    virtual void prefilter(TPrefilterCallback callback) override
    {
        callback(_databaseSegments, _prefilter_ptr->_queryFilter);
    }

private:
    NoQueryPrefilter * _prefilter_ptr;
    TDatabaseSegments _databaseSegments{};
};


} // namespace stellar
