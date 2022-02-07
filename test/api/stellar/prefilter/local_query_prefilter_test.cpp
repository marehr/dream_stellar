#include <gtest/gtest.h>

#include <seqan3/std/span>

// DEBUG: remove me
#include <seqan3/core/debug_stream.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <stellar/stellar_types.hpp>
#include <stellar/seqan2to3_alphabet.hpp>
#include <stellar/prefilter/detail/reverse_complementable_sequence.hpp>
#include <stellar/prefilter/detail/reinterpret_as_dna4_sequence.hpp>
#include <stellar/prefilter/detail/print_minimiser.hpp>
#include <stellar/prefilter/detail/valik_local_prefilter.hpp>
#include <stellar/prefilter/detail/valik_local_query_prefilter_index.hpp>

namespace stellar
{

// Prefilter: filters out the queries for each database sequence.
// In this case it does not filter out any queries and only builds a single query filter
template <typename TAlphabet, typename TQueryFilter_>
struct LocalQueryPrefilter
{
    using TQueryFilter = TQueryFilter_;

    LocalQueryPrefilter() = delete;

    LocalQueryPrefilter(
        StringSet<String<TAlphabet> > const & databases,
        StringSet<String<TAlphabet> > const & queries,
        stellar::detail::valik_local_query_prefilter_index queryPrefilterIndex
    ) :
        _databases{databases},
        _queries{queries},
        _queryPrefilterIndex{std::move(queryPrefilterIndex)}
    {}

    using TQuerySegment = seqan::Segment<seqan::String<TAlphabet> const, seqan::InfixSegment>;
    using TQuerySegments = std::span<TQuerySegment const>;

    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    using TDatabaseSegments = std::span<TDatabaseSegment const>;

    class Agent;

    std::vector<Agent> agents(size_t const agentCount, size_t const minLength)
    {
        assert(_querySegments.size() == 0u); // calling agents twice is forbidden

        // _databaseSegments = _splitter.split(_databases, minLength);
        // TODO: use SPLITTER
        for (String<TAlphabet> const & query: _queries)
            if (length(query) >= minLength)
                _querySegments.emplace_back(query, 0u, length(query));

        std::vector<Agent> agents{};
        agents.reserve(agentCount);

        // TODO: use DISTRIBUTER
        TQuerySegments allQuerySegments{_querySegments.begin(), _querySegments.end()};
        agents.emplace_back(*this, std::move(allQuerySegments), minLength);
        // for (TDatabaseSegments & segments: DatabaseAgentDistributer::distribute(_databaseSegments, agentCount))
        //     agents.emplace_back(*this, std::move(segments));

        return agents;
    }

private:
    StringSet<String<TAlphabet> > const & _databases{};
    StringSet<String<TAlphabet> > const & _queries{};
    std::vector<TQuerySegment> _querySegments{};
    stellar::detail::valik_local_query_prefilter_index _queryPrefilterIndex{};
};

template <typename TAlphabet, typename TQueryFilter>
struct LocalQueryPrefilter<TAlphabet, TQueryFilter>::Agent
{
    friend LocalQueryPrefilter;

    Agent(LocalQueryPrefilter & prefilter, TQuerySegments querySegments, size_t const minLength)
        : _prefilter_ptr{&prefilter}, _querySegments{std::move(querySegments)}, _minLength{minLength}
    {}
public:

    template <typename functor_t>
    void prefilter(functor_t && functor)
    {
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

            stellar::detail::valik_local_prefilter(
                query, // TODO: pass querySegment directly
                _prefilter_ptr->_queryPrefilterIndex,
                _minLength,
                [&](size_t const bin_id) {
                    _databaseSegments.emplace_back(binIDToDatabaseSegment(bin_id));
                });

            if (_databaseSegments.size() > 0)
            {
                TQueryFilter _queryFilter{querySegment};
                functor(_databaseSegments, _queryFilter);
            }
        }
    }

private:
    LocalQueryPrefilter * _prefilter_ptr{nullptr};
    TQuerySegments _querySegments{};
    size_t _minLength{};
};

} // namespace stellar

auto build_lqpindex()
{
    stellar::detail::valik_local_query_prefilter_index lqpindex
    {
        // ibf parameters
        ._ibf {
            seqan3::bin_count{5u /*arguments->bins*/},
            seqan3::bin_size{4096u / 1u/*arguments->bits / arguments->parts*/},
            seqan3::hash_function_count{2u/*arguments->hash*/}
        },

        // hash parameters
        ._kmer_size{seqan3::ungapped{4u/*arguments->kmer_size*/}},
        ._window_size{seqan3::window_size{6u/*arguments->window_size*/}}
    };

    return lqpindex;
}

TEST(LocalQueryPrefilter, LocalQueryPrefilter)
{
    using TAlphabet = seqan::Dna5;
    using TQueryFilter = seqan::Segment<seqan::String<TAlphabet> const, seqan::InfixSegment>;

    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "TAGAC" "TAACCGCAGAACACG" "ACTCCTCTAC");
    appendValue(databases, "TCAC" "GCCAGTTTAATTTTA" "GGGCTCCAAGT");
    appendValue(databases, "TACGCATATCTGG" "TAACCGCAGAACACG" "AA");

    auto lqpindex = build_lqpindex();
    for (auto const & database : databases)
    {
        size_t const i = &database - &databases[0];
        lqpindex.insert_sequence(seqan3::bin_index{i}, database);
    }

    seqan::StringSet<seqan::String<TAlphabet> > queries{};
    appendValue(queries, "TAGCCG" "TAACCGCAGAACACG" "GACC"); // matches database[0,2]
    appendValue(queries, "TA" "GCCAGTTTAATTTTA" "TACCAAGA"); // matches database[1]
    appendValue(queries, "ACAGCTACCGT"); // matches nothing (does not match minLength)
    appendValue(queries, "TTGAGACGAAT" "AGAACACGACTCCTCTAC" "CTTGTGGCA"); // matches database[0]

    // Note: TQueryFilter are the queries themselve to check whith which queries the actual filter would have been
    //       constructed.
    using TPrefilter = stellar::LocalQueryPrefilter<TAlphabet, TQueryFilter>;
    using TPrefilterAgent = TPrefilter::Agent;

    size_t const agentCount = 1u;
    size_t const minLength = 15u;
    TPrefilter prefilter{databases, queries, std::move(lqpindex)};

    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    std::vector<TDatabaseSegment> seenDatabases{};
    std::vector<TDatabaseSegment> expectedSeenDatabases
    {
        TDatabaseSegment{databases[0], 0, seqan::length(databases[0])}, // seen by queries[0]
        TDatabaseSegment{databases[0], 0, seqan::length(databases[0])}, // seen by queries[3]
        TDatabaseSegment{databases[1], 0, seqan::length(databases[1])}, // seen by queries[1]
        TDatabaseSegment{databases[2], 0, seqan::length(databases[2])}, // seen by queries[0]
    };

    using TQuerySegment = seqan::Segment<seqan::String<TAlphabet> const, seqan::InfixSegment>;
    std::vector<TQuerySegment> seenQueries{};
    std::vector<TQuerySegment> expectedSeenQueries
    {
        TQuerySegment{queries[0], 0, seqan::length(queries[0])},
        TQuerySegment{queries[1], 0, seqan::length(queries[1])},
        // will been filtered out due to beeing to short
        // TQuerySegment{queries[2], 0, seqan::length(queries[2])},
        TQuerySegment{queries[3], 0, seqan::length(queries[3])},
    };

    for (TPrefilterAgent prefilterAgent: prefilter.agents(agentCount, minLength))
    {
        using TDatabaseSegments = typename TPrefilter::TDatabaseSegments;
        using TQueryFilter = typename TPrefilter::TQueryFilter;

        prefilterAgent.prefilter([&](TDatabaseSegments databaseSegments, TQueryFilter queryFilter)
        {
#if STELLAR_DEBUG_MINIMISER
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << "queryID: " << (&host(queryFilter) - &queries[0]) << std::endl;
            std::cout << "querySegment: " << queryFilter << std::endl;
#endif // STELLAR_DEBUG_MINIMISER

            // collect all seen databases
            for (TDatabaseSegment databaseSegment : databaseSegments)
            {
#if STELLAR_DEBUG_MINIMISER
                std::cout << "\tdatabaseID: " << (&databaseSegment.underlyingDatabase() - &databases[0]) << std::endl;
                std::cout << "\tdatabaseSegment: " << databaseSegment.asInfixSegment() << std::endl;
#endif // STELLAR_DEBUG_MINIMISER
                seenDatabases.push_back(std::move(databaseSegment));
            }

            // collect all seen queries
            TQuerySegment querySegment = queryFilter; // NOTE: we choose queryFilter to be TQuerySegment
            seenQueries.push_back(std::move(querySegment));

#if STELLAR_DEBUG_MINIMISER
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            std::cout << std::endl;
#endif // STELLAR_DEBUG_MINIMISER
        });
    }

    // Note: every database should have been seen completely
    std::sort(seenDatabases.begin(), seenDatabases.end());
    std::sort(expectedSeenDatabases.begin(), expectedSeenDatabases.end());
    EXPECT_EQ(seenDatabases, expectedSeenDatabases);

    // Note: every query should have been seen completely
    std::sort(seenQueries.begin(), seenQueries.end());
    std::sort(expectedSeenQueries.begin(), expectedSeenQueries.end());
    EXPECT_EQ(seenQueries, expectedSeenQueries);
}
