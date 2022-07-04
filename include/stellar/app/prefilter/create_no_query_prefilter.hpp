
#pragma once

#include <stellar/stellar_types.hpp>

#include <stellar/app/prefilter/create_prefilter.hpp>

#include <stellar/prefilter/no_query_prefilter.hpp>
#include <stellar/prefilter/nsegment_database_agent_splitter.hpp>
#include <stellar/prefilter/whole_database_agent_splitter.hpp>

namespace stellar
{

namespace app
{

template <typename TAlphabet>
struct CreatePrefilter<NoQueryPrefilter<TAlphabet, NSegmentDatabaseAgentSplitter>>
{
    static std::unique_ptr<stellar::prefilter<TAlphabet>>
    create(
        StellarOptions const & options,
        StringSet<String<TAlphabet> > const & databases,
        StringSet<String<TAlphabet> > const & queries,
        StellarSwiftPattern<TAlphabet> & swiftPattern)
    {
        using TSplitter = NSegmentDatabaseAgentSplitter;
        using TPrefilter = NoQueryPrefilter<TAlphabet, TSplitter>;

        TSplitter splitter{options.splitNSegments};
        return std::make_unique<TPrefilter>(databases, swiftPattern /*copy pattern*/, splitter);
    }
};

template <typename TAlphabet>
struct CreatePrefilter<NoQueryPrefilter<TAlphabet, WholeDatabaseAgentSplitter>>
{
    static std::unique_ptr<stellar::prefilter<TAlphabet>>
    create(
        StellarOptions const & options,
        StringSet<String<TAlphabet> > const & databases,
        StringSet<String<TAlphabet> > const & queries,
        StellarSwiftPattern<TAlphabet> & swiftPattern)
    {
        using TSplitter = NSegmentDatabaseAgentSplitter;
        using TPrefilter = NoQueryPrefilter<TAlphabet, TSplitter>;

        TSplitter splitter{};
        return std::make_unique<TPrefilter>(databases, swiftPattern /*copy pattern*/, splitter);
    }
};

} // namespace app

} // namespace stellar
