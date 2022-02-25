
#pragma once

#include <stellar/stellar_types.hpp>

namespace stellar {

template <typename TVerifierTag>
struct SwiftHitVerifier
{
    using TEpsilon = double; // decltype(StellarOptions{}.epsilon)
    using TSize = int; // decltype(StellarOptions{}.minLength)
    using TDrop = double; // decltype(StellarOptions{}.xDrop)

    TEpsilon const epsilon;
    TSize const minLength;
    TDrop const xDrop;

    template <typename TAlphabet, typename TDelta, typename TOnAlignmentResultFn>
    void verify(StellarDatabaseSegment<TAlphabet> const & databaseSegment,
                StellarQuerySegment<TAlphabet> const & querySegment,
                TDelta const delta,
                TOnAlignmentResultFn && onAlignmentResult)
    {
        static_assert(std::is_signed<TSize>::value, "TDelta must be signed integral.");
        static_assert(std::is_floating_point<TDrop>::value, "TDrop must be floating point.");
        static_assert(std::is_unsigned<TDelta>::value, "TDelta must be unsigned integral.");

        this->verify(databaseSegment.asFinderSegment(), querySegment.asPatternSegment(), delta, onAlignmentResult);
    }

    template <typename TSequence, typename TDelta, typename TOnAlignmentResultFn>
    void verify(Segment<Segment<TSequence const, InfixSegment>, InfixSegment> const & finderSegment,
                Segment<Segment<TSequence const, InfixSegment>, InfixSegment> const & patternSegment,
                TDelta const delta,
                TOnAlignmentResultFn && onAlignmentResult)
    {
        static_assert(std::is_signed<TSize>::value, "TDelta must be signed integral.");
        static_assert(std::is_floating_point<TDrop>::value, "TDrop must be floating point.");
        static_assert(std::is_unsigned<TDelta>::value, "TDelta must be unsigned integral.");
        verifySwiftHit(finderSegment, patternSegment, epsilon, minLength, xDrop,
                       delta, onAlignmentResult, TVerifierTag{});
    }
};

} // namespace stellar
