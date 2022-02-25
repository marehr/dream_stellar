
#include <stellar/stellar_types.hpp>

#pragma once

namespace stellar::test
{

struct error_rate_t
{
    size_t alignment_length{};
    size_t error_columns{};
    double error_rate() const
    {
        return static_cast<double>(error_columns) / static_cast<double>(alignment_length);
    }
};

template <typename TSequence>
error_rate_t error_rate(seqan::Align<TSequence const, seqan::ArrayGaps> const & alignment)
{
    auto const & gappedSequence1 = row(alignment, 0);
    auto const & gappedSequence2 = row(alignment, 1);
    assert(length(gappedSequence1) == length(gappedSequence2));

    size_t alignment_length = length(gappedSequence1);
    size_t error_columns = 0u;
    for (size_t i = 0; i < alignment_length; ++i)
    {
        if (gappedSequence1[i] != gappedSequence2[i])
            ++error_columns;
    }

    error_rate_t error_rate
    {
        .alignment_length = alignment_length,
        .error_columns = error_columns,
    };
    return error_rate;
}

} // namespace stellar::test
