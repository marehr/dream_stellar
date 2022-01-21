
#pragma once

#include <seqan3/std/span>
#include <vector>

#include <stellar/stellar_database_segment.hpp>

namespace stellar
{

struct DatabaseAgentDistributer
{
    template <typename TAlphabet>
    using TStorage = std::vector<stellar::StellarDatabaseSegment<TAlphabet>>;

    template <typename TAlphabet>
    using TStorageSpan = std::span<stellar::StellarDatabaseSegment<TAlphabet> const>;

    template <typename TAlphabet>
    using TDistribution = std::vector<TStorageSpan<TAlphabet>>;

    template <typename TAlphabet>
    static inline TDistribution<TAlphabet> distribute(TStorage<TAlphabet> const & storage, size_t const threads)
    {
        TDistribution<TAlphabet> segmentsCollection{};
        segmentsCollection.reserve(threads + 1u);

        // split `storage` into `threads` many chunks
        size_t const segmentsPerThread = storage.size() / threads;
        size_t remainder = storage.size() % threads;
        size_t segmentLength = segmentsPerThread;

        for (size_t begin = 0u; begin < storage.size(); begin += segmentLength)
        {
            segmentLength = segmentsPerThread + (remainder > 0 ? (--remainder, 1) : 0);

            size_t const end = std::min(begin + segmentLength, storage.size());

            TStorageSpan<TAlphabet> segments{&storage[begin], &storage[end]};
            segmentsCollection.push_back(segments);
        }

        return segmentsCollection;
    }
};

} // namespace stellar
