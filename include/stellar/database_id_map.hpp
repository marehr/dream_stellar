
#pragma once

#include <stellar/stellar_types.hpp>

namespace stellar
{

template <typename TAlphabet, typename TId = CharString>
struct DatabaseIDMap
{
    size_t recordID(String<TAlphabet> const & database) const
    {
        String<TAlphabet> const * begin = &databases[0];
        String<TAlphabet> const * current = std::addressof(database);
        return current - begin;
    }

    TId const & databaseID(String<TAlphabet> const & database) const
    {
        return databaseIDs[recordID(database)];
    }

    StringSet<String<TAlphabet> > const & databases;
    StringSet<TId> const & databaseIDs;
};

} // namespace stellar
