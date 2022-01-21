
#pragma once

#include <stellar/stellar_types.hpp>

namespace stellar::test
{

template <typename TAlphabet>
bool expectStorage(std::vector<stellar::StellarDatabaseSegment<TAlphabet>> storage, size_t const minLength)
{
    std::sort(storage.begin(), storage.end());

    bool noErrors = true;

    // valid interval lengths
    for (size_t i = 0; i < storage.size(); ++i)
    {
        noErrors = noErrors && (storage[i].size() != 0);
        EXPECT_GT(storage[i].size(), 0u) << "storage[" << i << "]";

        noErrors = noErrors && (minLength <= storage[i].size());
        EXPECT_GE(storage[i].size(), minLength) << "storage[" << i << "]";

        noErrors = noErrors && (storage[i].size() <= length(storage[i].underlyingDatabase()));
        EXPECT_LE(storage[i].size(), length(storage[i].underlyingDatabase())) << "storage[" << i << "]";
    }

    // unique segments
    {
        std::vector<stellar::StellarDatabaseSegment<TAlphabet>> unique{storage};
        auto newLast = std::unique(unique.begin(), unique.end());

        noErrors = noErrors && (unique.end() == newLast); // no duplicate value
        EXPECT_EQ(unique.end(), newLast);
    }

    // do all segments cover the complete sequence?
    {
        using TKey = seqan::String<TAlphabet> const *;
        using TInterval = std::pair<std::size_t, std::size_t>;
        std::unordered_map<TKey, TInterval> intervalMap{};

        for (size_t i = 0; i < storage.size(); ++i)
        {
            TKey key = std::addressof(storage[i].underlyingDatabase());
            TInterval currentInterval = storage[i].interval();
            using TKeyValuePair = decltype(intervalMap)::value_type;

            auto intervalIt = intervalMap.find(key);
            if (intervalIt == intervalMap.end()) // key not contained
            {
                intervalMap.insert(TKeyValuePair{key, currentInterval});
            } else
            {
                TInterval & mapInterval = intervalMap.at(key);
                mapInterval.first = std::min(mapInterval.first, currentInterval.first);
                mapInterval.second = std::max(mapInterval.second, currentInterval.second);
            }
        }

        for (auto && [key, interval] : intervalMap)
        {
            seqan::String<TAlphabet> const & database = *key;

            noErrors = noErrors && (interval.first == 0u);
            EXPECT_EQ(interval.first, 0u);

            noErrors = noErrors && (interval.second == length(database));
            EXPECT_EQ(interval.second, length(database));
        }
    }

    return noErrors;
}

} // namespace stellar::test
