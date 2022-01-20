#pragma once

#include <stellar/stellar_database_segment.hpp>

namespace stellar
{

// keeps a single database as it is and does not split it any finer.
// multiple agents get work over a single
struct WholeDatabaseAgentSplitter
{
    template <typename TAlphabet>
    using TStorage = std::vector<stellar::StellarDatabaseSegment<TAlphabet>>;

    template <typename TAlphabet>
    TStorage<TAlphabet> split(seqan::StringSet<seqan::String<TAlphabet> > const & databases, size_t const minLength) const
    {
        TStorage<TAlphabet> storage{};
        storage.reserve(length(databases));

        // convert databases into database segments
        for (seqan::String<TAlphabet> const & database : databases)
            if (length(database) >= minLength)
                storage.emplace_back(database, 0u, length(database));

        return storage;
    }
};

} // namespace stellar
