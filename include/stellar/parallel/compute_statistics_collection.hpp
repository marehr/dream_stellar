
#pragma once

#include <stellar/stellar_types.hpp>

namespace stellar
{

struct StellarComputeStatisticsCollection;

struct StellarComputeStatisticsPartialCollection
{
    friend StellarComputeStatisticsCollection;

    explicit StellarComputeStatisticsPartialCollection(size_t const databasesCount)
        : _statistics(databasesCount)
    {}

    void updateByRecordID(size_t const databaseRecordID, StellarComputeStatistics const & statistics)
    {
        _statistics[databaseRecordID].mergeIn(statistics);
    }
private:
    std::vector<StellarComputeStatistics> _statistics; // one per database
};

struct StellarComputeStatisticsCollection
{
    explicit StellarComputeStatisticsCollection(size_t const databasesCount)
        : _statistics(databasesCount)
    {}

    StellarComputeStatistics const & operator[](size_t const databaseRecordID) const
    {
        return _statistics[databaseRecordID];
    }

    void mergePartialIn(StellarComputeStatisticsPartialCollection const & partialStatistics)
    {
        for (size_t i = 0; i < _statistics.size(); ++i)
            _statistics[i].mergeIn(partialStatistics._statistics[i]);
    }

    size_t size() const
    {
        return _statistics.size();
    }
private:
    std::vector<StellarComputeStatistics> _statistics; // one per database
};

} // namespace stellar
