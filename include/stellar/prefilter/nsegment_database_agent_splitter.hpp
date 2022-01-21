#pragma once

#include <seqan3/std/span>

#include <stellar/stellar_database_segment.hpp>

namespace stellar
{

// keeps a single database as it is and does not split it any finer.
// multiple agents get work over a single
struct NSegmentDatabaseAgentSplitter
{
    template <typename TAlphabet>
    using TStorage = std::vector<stellar::StellarDatabaseSegment<TAlphabet>>;

    template <typename TAlphabet>
    TStorage<TAlphabet> split(seqan::StringSet<seqan::String<TAlphabet> > const & databases, size_t const minLength) const
    {
        assert(segmentCount > 0u);

        TStorage<TAlphabet> storage{};
        storage.reserve(length(databases));

        // convert databases into database segments
        for (seqan::String<TAlphabet> const & database : databases)
        {
            size_t const maxLength = length(database);
            size_t const maxSegmentCount = (minLength > 0u) ? maxLength / minLength : ~size_t{0u};
            size_t const segmentCount = std::min(this->segmentCount, maxSegmentCount);

            // segmentCount means maxLength < minLength
            if (segmentCount == 0)
                continue;

            size_t const minSegmentLength = maxLength / segmentCount;
            size_t segmentLength = minSegmentLength;
            size_t remainder = maxLength - segmentLength * segmentCount;

            for (size_t begin = 0u; begin < maxLength; begin += segmentLength)
            {
                segmentLength = minSegmentLength + (remainder > 0 ? (--remainder, 1) : 0);
                size_t const end = std::min(begin + segmentLength, maxLength);

                if (begin > 0u && minLength > 0u)
                    storage.emplace_back(database, begin - minLength, begin + minLength);
                storage.emplace_back(database, begin, end);
            }
        }

        return storage;
    }

    size_t segmentCount{1u};
};

} // namespace stellar
