
#pragma once

#include <stellar/stellar_sequence_segment.hpp>

namespace stellar
{

template <typename TAlphabet>
struct StellarDatabaseSegment : public StellarSequenceSegment<TAlphabet>
{
    using TBase = StellarSequenceSegment<TAlphabet>;

    using typename TBase::TInfixSegment;
    using TNestedFinderSegment = seqan::Segment<TInfixSegment, seqan::InfixSegment>;

    using TBase::TBase; // import constructor

    static StellarDatabaseSegment<TAlphabet> fromFinderMatch(TInfixSegment const & finderMatch)
    {
        seqan::String<TAlphabet> const & underlyingDatabase = host(finderMatch);
        return {underlyingDatabase, seqan::beginPosition(finderMatch), seqan::endPosition(finderMatch)};
    }

    seqan::String<TAlphabet> const & underlyingDatabase() const &
    {
        return this->underlyingSequence();
    }

    TNestedFinderSegment asFinderSegment() const
    {
        seqan::String<TAlphabet> const & _database = underlyingDatabase();
        auto finderInfix = this->asInfixSegment();

        TInfixSegment const finderInfixSeq = infix(_database, 0, length(_database));
        TNestedFinderSegment finderSegment(finderInfixSeq,
            seqan::beginPosition(finderInfix) - seqan::beginPosition(_database),
            seqan::endPosition(finderInfix) - seqan::beginPosition(_database));
        return finderSegment;
    }
};

} // namespace stellar
