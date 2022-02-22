// ==========================================================================
//                    STELLAR - SwifT Exact LocaL AligneR
//                   http://www.seqan.de/projects/stellar/
// ==========================================================================
// Copyright (C) 2010-2012 by Birte Kehr
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your options) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STELLAR_H
#define SEQAN_HEADER_STELLAR_H

#include <iostream>
#include <seqan/seeds.h>

#include <stellar/stellar_types.hpp>
#include <stellar/stellar_extension.hpp>
#include <stellar/stellar_index.hpp>
#include <stellar/verification/all_local.hpp>
#include <stellar/verification/banded_global_extend.hpp>
#include <stellar/verification/banded_global.hpp>
#include <stellar/verification/best_local.hpp>
#include <stellar/verification/swift_hit_verifier.hpp>

#include <stellar/app/stellar.diagnostics.hpp>

namespace stellar
{

using namespace seqan;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Checks whether two matches overlap *in seq2* and
//  whether the non-overlaping parts are shorter than minLength.
template<typename TMatch, typename TSize>
bool
checkOverlap(TMatch const & matchA, TMatch const & matchB, TSize const minLength) {
    // check id and orienation
    if (matchA.id != matchB.id || matchA.orientation != matchB.orientation) return false;
    if (matchA.id == TMatch::INVALID_ID || matchB.id == TMatch::INVALID_ID) return false;

    // check overlap in seq2
    if (matchA.begin2 >= matchB.begin2) {
        if (matchA.end2 >= matchB.end2) {
            // check length of non-overlapping parts of both matches
            if ((TSize)matchA.begin2 - (TSize)matchB.begin2 >= minLength &&
                (TSize)matchA.end2 - (TSize)matchB.end2 >= minLength) {
                return false;
            }
        }
        // check whether offset is the same in both sequences
        if (toViewPosition(matchA.row2, matchA.begin2) - toViewPosition(matchB.row2, matchB.begin2) !=
            toViewPosition(matchA.row1, matchA.begin1) - toViewPosition(matchB.row1, matchB.begin1)) {
            return false;
        }
    } else {
        if (matchA.end2 < matchB.end2) {
            // check length of non-overlapping parts of both matches
            if ((TSize)matchB.begin2 - (TSize)matchA.begin2 >= minLength &&
                (TSize)matchB.end2 - (TSize)matchA.end2 >= minLength) {
                return false;
            }
        }
        // check whether offset is the same in both sequences
        if (toViewPosition(matchB.row2, matchB.begin2) - toViewPosition(matchA.row2, matchA.begin2) !=
            toViewPosition(matchB.row1, matchB.begin1) - toViewPosition(matchA.row1, matchA.begin1)) {
            return false;
        }
    }
    return true;
}

template<typename TRow, typename TPosition>
TPosition
projectedPosition(TRow const & rowA, TRow const & rowB, TPosition pos)
{
    return toSourcePosition(rowB, toViewPosition(rowA, pos));
}

///////////////////////////////////////////////////////////////////////////////
// Checks all alignment columns of two overlapping matches.
// It is assumed that matchA.begin1 < matchB.begin1.
template<typename TMatch, typename TSize>
bool
_checkAlignColOverlap(TMatch const & matchA, TMatch const & matchB, TSize const minLength)
{
    TSize equalCols = 0;
    TSize diffCols = 0;

    for  (typename TMatch::TPos pos = matchB.begin1; pos < _min(matchA.end1, matchB.end1); ++pos)
    {
        if (projectedPosition(matchA.row1, matchA.row2, pos) == projectedPosition(matchB.row1, matchB.row2, pos))
            ++equalCols;
        else
            ++diffCols;
    }

    if (diffCols >= minLength) return false;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
// Marks matches that overlap in both sequences with a longer match as invalid.
template<typename TSequence, typename TId, typename TSize>
void maskOverlaps(String<StellarMatch<TSequence const, TId> > & matches, TSize const minLength)
{
    typedef StellarMatch<TSequence const, TId>              TMatch;
    typedef typename TMatch::TPos                           TPos;
    typedef typename Iterator<String<TMatch>, Rooted>::Type TIter;
    typedef typename Iterator<String<TSize>, Rooted>::Type  TOverlapIter;

    // sort matches by begin position in row0
    sortMatches(matches, LessPos<TMatch>());

    // iterate all matches
    TIter it = begin(matches);

    // list of matches that potentially overlap with the current match in row0 and
    // start earlier (including matches that overlap but have a unique part of at
    // least minLength) sorted by descending end positions
    String<TSize> overlaps;

    for (; it != end(matches); ++it)
    {
        if ((*it).id == TMatch::INVALID_ID) continue;

        TPos insertPos = 0;

        // iterate potentially overlapping matches
        TOverlapIter overlapIt = begin(overlaps);
        for (; overlapIt != end(overlaps); ++overlapIt)
        {
            TMatch & o = matches[*overlapIt];

            // determine position for inserting *it into overlaps after checking
            if ((*it).end1 < o.end1) insertPos++;

            // check if matches overlap in row0 - if not, then break
            if (o.end1 <= (*it).begin1) break;

            // check if unique parts of the two matches in row0 are longer than minLength - if yes, then continue
            if ((*it).begin1 - o.begin1 >= (TPos)minLength &&
                (*it).end1 > o.end1 && (*it).end1 - o.end1 >= (TPos)minLength) continue;

            // check if matches overlap in row1 - if not, then continue
            if (!checkOverlap(*it, o, minLength)) continue;

            // check exact alignment columns for overlap
            if (!_checkAlignColOverlap(o, *it, minLength)) continue;

            // set shorter match invalid
            if (length(*it) > length(o))
                o.id = TMatch::INVALID_ID;
            else
                (*it).id = TMatch::INVALID_ID;
        }

        // remove all matches from overlaps that end earlier than current match begins
        resize(overlaps, position(overlapIt));

        if ((*it).id != TMatch::INVALID_ID)
            insertValue(overlaps, insertPos, position(it));
    }

}

///////////////////////////////////////////////////////////////////////////////
// Removes matches that are marked as invalid, and then keeps only the numMatches best matches.
template<typename TSequence, typename TId, typename TSize>
void
compactMatches(String<StellarMatch<TSequence const, TId> > & matches, TSize const numMatches) {
    typedef StellarMatch<TSequence const, TId>                  TMatch;
    typedef typename Iterator<String<TMatch>, Standard>::Type   TIterator;

    // sort matches by length (and validity)
    sortMatches(matches, LessLength<TMatch>());

    // count valid matches
    TSize num = 0;

    TIterator it = begin(matches, Standard());
    TIterator itEnd = end(matches, Standard());

    for(; it != itEnd; ++it) {
        if ((*it).id != TMatch::INVALID_ID)
            ++num;
    }

    // keep only valid and longest matches
    resize(matches, _min(num, numMatches));
}

template<typename TMatch_>
inline bool
QueryMatches<TMatch_>::
removeOverlapsAndCompactMatches(size_t const disableThresh,
                                size_t const compactThresh,
                                size_t const minLength,
                                size_t const numMatches)
{
    if (this->disabled)
        return false;

    size_t const matchesCount = length(this->matches);

    if (matchesCount > disableThresh) {
        this->disabled = true;
        clear(this->matches);
        return false;
    }

    if (matchesCount <= compactThresh)
        return false;

    maskOverlaps(this->matches, minLength);      // remove overlaps and duplicates
    compactMatches(this->matches, numMatches);   // keep only the <numMatches> longest matches
    return true;
}

///////////////////////////////////////////////////////////////////////////////
// Appends a match to matches container and removes overlapping matches if threshold is reached.
template<typename TSource, typename TId, typename TSize, typename TSize1>
inline bool
_insertMatch(QueryMatches<StellarMatch<TSource const, TId> > & queryMatches,
             StellarMatch<TSource const, TId> const & match,
             TSize const minLength,
             TSize1 const disableThresh,
             TSize1 & compactThresh,
             TSize1 const numMatches) {

    appendValue(queryMatches.matches, match);

    // std::cerr << "Inserting match \n-------------\n" << match.row1 <<"\n" << match.row2 << "----------------\n";

    if (queryMatches.removeOverlapsAndCompactMatches(disableThresh, compactThresh, minLength, numMatches))
    {
        // raise compact threshold if many matches are kept
        if ((length(queryMatches.matches) << 1) > compactThresh)
            compactThresh += (compactThresh >> 1);
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////////
// Calls swift filter and verifies swift hits. = Computes eps-matches.
// A basic block for stellar
template<typename TAlphabet, typename TSource, typename TId, typename TTag>
StellarComputeStatistics
_stellarKernel(StellarSwiftFinder<TAlphabet> & finder,
               StellarSwiftPattern<TAlphabet> & pattern,
               StringSet<QueryMatches<StellarMatch<TSource const, TId> > > & matches,
               SwiftHitVerifier<TTag> & swiftVerifier) {
    using TText = String<TAlphabet>;
    using TMatch = StellarMatch<TSource const, TId>;
    using TInfix = typename Infix<TText const>::Type;

    StellarComputeStatistics statistics{};

    while (find(finder, pattern, swiftVerifier.epsilon, swiftVerifier.minLength)) {
        TInfix const finderInfix = infix(finder);
        TInfix const finderInfixSeq = infix(haystack(finder), 0, length(haystack(finder)));
        Segment<TInfix, InfixSegment> finderSegment(finderInfixSeq,
            beginPosition(finderInfix) - beginPosition(haystack(finder)),
            endPosition(finderInfix) - beginPosition(haystack(finder)));

        ++statistics.numSwiftHits;
        statistics.totalLength += length(finderInfix);
        statistics.maxLength = std::max<size_t>(statistics.maxLength, length(finderInfix));

        QueryMatches<TMatch> & queryMatches = value(matches, pattern.curSeqNo);

        if (queryMatches.disabled) continue;

        TText const & patternSeq = getSequenceByNo(pattern.curSeqNo, indexText(needle(pattern)));
        TInfix const patternInfix = infix(pattern, patternSeq);
        TInfix const patternInfixSeq = infix(patternSeq, 0, length(patternSeq));
        Segment<TInfix, InfixSegment> patternSegment(patternInfixSeq,
            beginPosition(patternInfix) - beginPosition(patternSeq),
            endPosition(patternInfix) - beginPosition(patternSeq));

        ////Debug stuff:
        //std::cout << beginPosition(finderInfix) << ",";
        //std::cout << endPosition(finderInfix) << "  ";
        //std::cout << beginPosition(patternSegment) << ",";
        //std::cout << endPosition(patternSegment) << std::endl;

        // verification
        swiftVerifier.verify(finderSegment,
                             patternSegment,
                             pattern.bucketParams[0].delta + pattern.bucketParams[0].overlap,
                             queryMatches);
    }

    return statistics;
}

///////////////////////////////////////////////////////////////////////////////
// Calls swift filter and verifies swift hits. = Computes eps-matches.
template<typename TAlphabet, typename TSize, typename TDrop, typename TSize1,
         typename TMode, typename TSource, typename TId, typename TTag>
void stellar(StellarSwiftFinder<TAlphabet> & finder,
             StellarSwiftPattern<TAlphabet> & pattern,
             double const epsilon,
             TSize const minLength,
             TDrop const xDrop,
             TSize1 const disableThresh,
             TSize1 & compactThresh,
             TSize1 const numMatches,
             TMode const verbose,
             TId const & databaseID,
             bool const databaseStrand,
             StringSet<QueryMatches<StellarMatch<TSource const, TId> > > & matches,
             TTag /*tag*/) {

    SwiftHitVerifier<TTag> swiftVerifier
    {
        STELLAR_DESIGNATED_INITIALIZER(.epsilon =, epsilon),
        STELLAR_DESIGNATED_INITIALIZER(.minLength =, minLength),
        STELLAR_DESIGNATED_INITIALIZER(.xDrop =, xDrop),
        STELLAR_DESIGNATED_INITIALIZER(.disableThresh =, disableThresh),
        STELLAR_DESIGNATED_INITIALIZER(.compactThresh =, compactThresh),
        STELLAR_DESIGNATED_INITIALIZER(.numMatches =, numMatches),
        STELLAR_DESIGNATED_INITIALIZER(.databaseID =, databaseID),
        STELLAR_DESIGNATED_INITIALIZER(.databaseStrand =, databaseStrand)
    };

    StellarComputeStatistics statistics = _stellarKernel(finder, pattern, matches, swiftVerifier);

    if (verbose)
        stellar::app::_printStellarKernelStatistics(statistics);

    using TMatch = StellarMatch<TSource const, TId>;
    for (QueryMatches<TMatch> & queryMatches : matches) {
        removeOverlapsAndCompactMatches(queryMatches, disableThresh, /*compactThresh*/ 0, minLength, numMatches);
    }
}

// Wrapper for stellar
template<typename TAlphabet, typename TSize, typename TDrop,
         typename TSource, typename TId, typename TTag>
void stellar(StellarSwiftFinder<TAlphabet> & finder,
             StellarSwiftPattern<TAlphabet> & pattern,
             double const epsilon,
             TSize const minLength,
             TDrop const xDrop,
             StringSet<QueryMatches<StellarMatch<TSource const, TId> > > & matches,
             TTag tag) {
    unsigned const maxValue = (unsigned)-1;
    unsigned const disableThresh = maxValue;
    unsigned const numMatches = disableThresh;
    unsigned compactThresh = 1000;
    TId const id = "db";

    stellar(finder, pattern, epsilon, minLength, xDrop,
            disableThresh, compactThresh, numMatches, /*verbose*/ false, id, /*databaseStrand*/ true,
            matches, tag);
}

} // namespace stellar

#endif