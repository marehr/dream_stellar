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

#ifndef SEQAN_HEADER_STELLAR_TYPES_H
#define SEQAN_HEADER_STELLAR_TYPES_H

#include <seqan/align.h>

#if __cpp_designated_initializers || __GNUC__ >= 8
#   define STELLAR_DESIGNATED_INITIALIZER(designator, value) designator value
#else
#   define STELLAR_DESIGNATED_INITIALIZER(designator, value) value
#endif // __cpp_designated_initializers || __GNUC__ >= 8

namespace stellar
{

using namespace seqan;

struct VerifyAllLocal_;
typedef Tag<VerifyAllLocal_> const AllLocal;

struct VerifyBestLocal_;
typedef Tag<VerifyBestLocal_> const BestLocal;

struct VerifyBandedGlobal_;
typedef Tag<VerifyBandedGlobal_> const BandedGlobal;

struct VerifyBandedGlobalExtend_;
typedef Tag<VerifyBandedGlobalExtend_> const BandedGlobalExtend;

// basically a std::variant<AllLocal, BestLocal, BandedGlobal, BandedGlobalExtend>
struct StellarVerificationMethod
{
    StellarVerificationMethod(AllLocal) : _index{0} {}
    StellarVerificationMethod(BestLocal) : _index{1} {}
    StellarVerificationMethod(BandedGlobal) : _index{2} {}
    StellarVerificationMethod(BandedGlobalExtend) : _index{3} {}

    constexpr std::size_t index() const noexcept
    {
        return _index;
    }

    friend constexpr bool operator==(
        StellarVerificationMethod const & m1,
        StellarVerificationMethod const & m2)
    {
        return m1.index() == m2.index();
    }

    friend inline std::string to_string(StellarVerificationMethod method)
    {
        using cstring_t = char const * const;
        cstring_t method_names[] = {"exact", "bestLocal", "bandedGlobal", "bandedGlobalExtend"};
        return method_names[method.index()];
    }

private:
    std::size_t _index{~std::size_t{0u}};
};

///////////////////////////////////////////////////////////////////////////////
// Options for Stellar
struct StellarOptions {
    // i/o options
    CharString databaseFile;        // name of database file
    CharString queryFile;           // name of query file
    CharString outputFile;          // name of result file
    CharString disabledQueriesFile; // name of result file containing disabled queries
    CharString outputFormat;        // Possible formats: gff, text
    CharString alphabet;            // Possible values: dna, rna, protein, char
    bool noRT;                      // suppress printing of running time if set to true

    // main options
    unsigned qGram;             // length of the q-grams
    double epsilon;             // maximal error rate
    int minLength;              // minimal length of an epsilon-match
    double xDrop;               // maximal x-drop

    // more options
    unsigned threadCount{1u};   // The maximum number of threads
    bool splitDatabase{true};   // Whether the database sequences should be split into multiple parts.
                                // This can result into a slight increase in the number of Swift verifications.
    bool forward;               // compute matches to forward strand of database
    bool reverse;               // compute matches to reverse complemented database
    // verification strategy: exact, bestLocal, bandedGlobal
    StellarVerificationMethod verificationMethod{AllLocal{}};
    unsigned disableThresh;     // maximal number of matches allowed per query before disabling verification of hits for that query
    unsigned compactThresh;     // number of matches after which removal of overlaps and duplicates is started
    unsigned numMatches;        // maximal number of matches per query and database
    unsigned maxRepeatPeriod;   // maximal period of low complexity repeats to be filtered
    unsigned minRepeatLength;   // minimal length of low complexity repeats to be filtered
    double qgramAbundanceCut;
    bool verbose;               // verbose mode


    StellarOptions() {
        outputFile = "stellar.gff";
        disabledQueriesFile = "stellar.disabled.fasta";
        outputFormat = "gff";
        alphabet = "dna5";
        noRT = false;

        qGram = std::numeric_limits<unsigned>::max();
        epsilon = 0.05;
        minLength = 100;
        xDrop = 5;

        forward = true;
        reverse = true;
        disableThresh = std::numeric_limits<unsigned>::max();
        compactThresh = 500;
        numMatches = 50;
        maxRepeatPeriod = 1;
        minRepeatLength = 1000;
        qgramAbundanceCut = 1;
        verbose = false;
    }
};

struct StellarStatistics
{
    bool kMerComputed{};
    unsigned kMerLength{};
    unsigned smin{};
    int threshold{};
    int distanceCut{};
    int delta{};
    int overlap{};

    StellarStatistics(StellarOptions const & options)
    {

        int errMinLen = (int) floor(options.epsilon * options.minLength);
        int n = (int) ceil((errMinLen + 1) / options.epsilon);
        int errN = (int) floor(options.epsilon * n);
        smin = (unsigned) _min(ceil((double)(options.minLength - errMinLen) / (errMinLen + 1)),
                               ceil((double)(n - errN) / (errN + 1)));

        kMerLength = options.qGram;
        kMerComputed = options.qGram == (unsigned)-1;

        if (kMerComputed)
            kMerLength = (unsigned)_min(smin, 32u);

        threshold = (int) _max(1, (int) _min((n + 1) - options.qGram * (errN + 1),
                                             (options.minLength + 1) - options.qGram * (errMinLen + 1)));
        overlap = (int) floor((2 * threshold + options.qGram - 3) / (1 / options.epsilon - options.qGram));
        distanceCut = (threshold - 1) + options.qGram * overlap + options.qGram;
        int logDelta = _max(4, (int) ceil(log((double)overlap + 1) / log(2.0)));
        delta = 1 << logDelta;
    }
};

struct StellarComputeStatistics
{
    size_t numSwiftHits = 0;

    size_t maxLength = 0;
    size_t totalLength = 0;

    void mergeIn(StellarComputeStatistics const & statistics)
    {
        this->numSwiftHits += statistics.numSwiftHits;
        this->totalLength += statistics.totalLength;
        this->maxLength = std::max<size_t>(this->maxLength, statistics.maxLength);
    }
};

struct StellarOutputStatistics
{
    size_t maxLength{0u};
    size_t totalLength{0u};
    size_t numMatches{0u};
    size_t numDisabled{0u};

    void mergeIn(StellarOutputStatistics const & statistics)
    {
        maxLength = std::max(maxLength, statistics.maxLength);
        totalLength = totalLength + statistics.totalLength;
        numMatches = numMatches + statistics.numMatches;
        numDisabled = numDisabled + statistics.numDisabled;
    }
};

///////////////////////////////////////////////////////////////////////////////
// Container for storing local alignment matches of one query sequence
template<typename TMatch_>
struct QueryMatches {
    typedef typename Size<typename Source<TMatch_>::Type>::Type TSize;

    String<TMatch_> matches;
    bool disabled;
    TSize lengthAdjustment;

    QueryMatches() : disabled(false), lengthAdjustment(0)
    {}

    bool removeOverlapsAndCompactMatches(size_t const disableThresh,
                                         size_t const compactThresh,
                                         size_t const minLength,
                                         size_t const numMatches);

    void mergeIn(QueryMatches const & otherMatches)
    {
        this->disabled = this->disabled || otherMatches.disabled;
        if (this->disabled)
        {
            clear(this->matches);
        } else
        {
            append(this->matches, otherMatches.matches);
        }
    }
};

///////////////////////////////////////////////////////////////////////////////
// Container for storing a local alignment match
template<typename TSequence_, typename TId_>
struct StellarMatch {
    static_assert(std::is_const<TSequence_>::value, "Sequence must be const qualified! I.e. StellarMatch<... const, ...>");
    typedef TSequence_                          TSequence;
    typedef TId_                                TId;
    typedef typename Position<TSequence>::Type  TPos;

    typedef Align<TSequence, ArrayGaps>         TAlign;
    typedef typename Row<TAlign>::Type         TRow;

    static const TId INVALID_ID;

    TId id;         // database ID
    bool orientation;
    TPos begin1;
    TPos end1;
    TRow row1;

    TPos begin2;
    TPos end2;
    TRow row2;

    StellarMatch() : id(), orientation(false), begin1(0), end1(0), begin2(0), end2(0)
    {}

    template <typename TAlign, typename TId>
    StellarMatch(TAlign & _align, TId _id, bool _orientation)
    {
        id = _id;
        orientation = _orientation;
        row1 = row(_align, 0);
        row2 = row(_align, 1);

        begin1 = beginPosition(row1);
        end1 = endPosition(row1);

        begin2 = beginPosition(row2);
        end2 = endPosition(row2);
    }
};

template <typename TSequence, typename TId>
const TId
StellarMatch<TSequence, TId>::INVALID_ID = "###########";


///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// to sort matches by position and remove overlapping matches
template <typename TMatch>
struct LessPos : public ::std::function<bool(TMatch, TMatch)> {
    LessPos() {}

    inline int compare(TMatch const & a, TMatch const & b) const {
        // query number
        if ((a.id) < (b.id)) return -1;
        if ((a.id) > (b.id)) return 1;

        // database begin position
        typename TMatch::TPos aBegin1 = _min(a.begin1, a.end1);
        typename TMatch::TPos bBegin1 = _min(b.begin1, b.end1);
        if (aBegin1 < bBegin1) return -1;
        if (aBegin1 > bBegin1) return 1;

        // database end position
        typename TMatch::TPos aEnd1 = _max(a.begin1, a.end1);
        typename TMatch::TPos bEnd1 = _max(b.begin1, b.end1);
        if (aEnd1 < bEnd1) return -1;
        if (aEnd1 > bEnd1) return 1;

        // query begin position
        typename TMatch::TPos aBegin2 = _min(a.begin2, a.end2);
        typename TMatch::TPos bBegin2 = _min(b.begin2, b.end2);
        if (aBegin2 < bBegin2) return -1;
        if (aBegin2 > bBegin2) return 1;

        // query end position
        typename TMatch::TPos aEnd2 = _max(a.begin2, a.end2);
        typename TMatch::TPos bEnd2 = _max(b.begin2, b.end2);
        if (aEnd2 < bEnd2) return -1;
        if (aEnd2 > bEnd2) return 1;

        // orientation
        if (a.orientation > b.orientation) return -1;
        if (a.orientation < b.orientation) return 1;

        //// orientation
        //bool oa = a.begin1 < a.end1;
        //bool ob = b.begin1 < b.end1;
        //if (oa != ob) return oa;

        return 0;
    }

    inline bool operator() (TMatch const & a, TMatch const & b) const {
        return compare(a, b) == -1;
    }
};

///////////////////////////////////////////////////////////////////////////////
// to sort matches by length
template <typename TMatch>
struct LessLength : public ::std::function<bool(TMatch, TMatch)> {
    LessLength() {}

    inline int compare(TMatch const & a, TMatch const & b) const {
        if (a.id == TMatch::INVALID_ID) return 1;
        if (b.id == TMatch::INVALID_ID) return -1;

        typename TMatch::TPos aLength = abs((int)a.end1 - (int)a.begin1);
        typename TMatch::TPos bLength = abs((int)b.end1 - (int)b.begin1);

        if (aLength < bLength) return 1;
        if (aLength > bLength) return -1;

        return 0;
    }

    inline bool operator() (TMatch const & a, TMatch const & b) const {
        return compare(a, b) == -1;
    }
};

///////////////////////////////////////////////////////////////////////////////
// Determines whether match1 is upstream of match2 in specified row.
//  If matches overlap, the non-overlapping parts have to be longer than minLenght.
template<typename TSequence, typename TId, typename TRowNo, typename TSize>
inline bool
_isUpstream(StellarMatch<TSequence const, TId> & match1, StellarMatch<TSequence const, TId> & match2, TRowNo row, TSize minLength) {
    typedef typename StellarMatch<TSequence const, TId>::TPos TPos;

    TPos e1, b2;
    if (row == 0) {
        e1 = match1.end1;
        b2 = match2.begin1;
    } else {
        e1 = match1.end2;
        b2 = match2.begin2;
    }

    if (e1 <= b2) return true;

    TPos b1, e2;
    if (row == 0) {
        e2 = match2.end1;
        b1 = match1.begin1;
    } else {
        e2 = match2.end2;
        b1 = match1.begin2;
    }

    if (b1 < b2 && (b2 - b1 >= minLength)) {
        if ((e1 < e2) && (e2 - e1 >= minLength)) return true;
    }

    return false;
}

///////////////////////////////////////////////////////////////////////////////
// sorts StellarMatchees by specified functor
template <typename TSequence, typename TId, typename TFunctorLess>
inline void
sortMatches(String<StellarMatch<TSequence const, TId> > & stellarMatches, TFunctorLess const & less) {
    std::stable_sort(
        begin(stellarMatches, Standard()),
        end(stellarMatches, Standard()),
        less);
}

///////////////////////////////////////////////////////////////////////////////
// returns the length of the longer row from StellarMatch
template <typename TSequence, typename TId>
inline typename Size<TSequence const>::Type
length(StellarMatch<TSequence const, TId> & match) {
    return _max(length(match.row1), length(match.row2));
}

} // namespace stellar

#endif
