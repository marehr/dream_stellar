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

#pragma once

#include <stellar/app/stellar.main.hpp>

#include <seqan/seq_io.h>

#include <stellar/stellar.hpp>
#include <stellar/stellar_index.hpp>
#include <stellar/stellar_output.hpp>

#include <stellar/parallel/compute_statistics_collection.hpp>

#include <stellar/app/stellar.diagnostics.hpp>

#include <stellar/app/stellar.diagnostics.tpp>

namespace stellar
{

using namespace seqan;

namespace app
{

///////////////////////////////////////////////////////////////////////////////
// Initializes a Finder object for a database sequence,
//  calls stellar, and writes matches to file
template <typename TAlphabet, typename TIsPatternDisabledFn, typename TOnAlignmentResultFn>
inline StellarComputeStatistics
_stellarOnOne(StellarDatabaseSegment<TAlphabet> const & databaseSegment,
              StellarSwiftPattern<TAlphabet> & swiftPattern,
              StellarOptions const & options,
              TIsPatternDisabledFn && isPatternDisabled,
              TOnAlignmentResultFn && onAlignmentResult)
{
    // finder
    StellarSwiftFinder<TAlphabet> swiftFinder(databaseSegment.asInfixSegment(), options.minRepeatLength, options.maxRepeatPeriod);

    auto _stellar = [&](auto tag) -> StellarComputeStatistics
    {
        using TTag = decltype(tag);
        SwiftHitVerifier<TTag> swiftVerifier
        {
            STELLAR_DESIGNATED_INITIALIZER(.epsilon = , options.epsilon),
            STELLAR_DESIGNATED_INITIALIZER(.minLength = , options.minLength),
            STELLAR_DESIGNATED_INITIALIZER(.xDrop = , options.xDrop)
        };

        return _stellarKernel(swiftFinder, swiftPattern, swiftVerifier, isPatternDisabled, onAlignmentResult);
    };

    StellarComputeStatistics statistics;

    // stellar
    if (options.verificationMethod == StellarVerificationMethod{AllLocal{}})
        statistics = _stellar(AllLocal());
    else if (options.verificationMethod == StellarVerificationMethod{BestLocal{}})
        statistics = _stellar(BestLocal());
    else if (options.verificationMethod == StellarVerificationMethod{BandedGlobal{}})
        statistics = _stellar(BandedGlobal());
    else if (options.verificationMethod == StellarVerificationMethod{BandedGlobalExtend{}})
        statistics = _stellar(BandedGlobalExtend());

    return statistics;
}

} // namespace stellar::app

} // namespace stellar

//////////////////////////////////////////////////////////////////////////////
namespace seqan {

template <typename TAlphabet>
struct Cargo<::stellar::StellarQGramIndex<TAlphabet>>
{
    typedef struct
    {
        double      abundanceCut;
    } Type;
};

//////////////////////////////////////////////////////////////////////////////
// Repeat masker
template <typename TAlphabet>
inline bool _qgramDisableBuckets(::stellar::StellarQGramIndex<TAlphabet> & index)
{
    using TIndex = ::stellar::StellarQGramIndex<TAlphabet>;
    using TDir = typename Fibre<TIndex, QGramDir>::Type;
    using TDirIterator = typename Iterator<TDir, Standard>::Type;
    using TSize = typename Value<TDir>::Type;

    TDir & dir   = indexDir(index);
    bool result  = false;
    unsigned counter = 0;
    TSize thresh = (TSize)(length(index) * cargo(index).abundanceCut);
    if (thresh < 100)
        thresh = 100;

    TDirIterator it = begin(dir, Standard());
    TDirIterator itEnd = end(dir, Standard());
    for (; it != itEnd; ++it)
        if (*it > thresh)
        {
            *it = (TSize) - 1;
            result = true;
            ++counter;
        }

    if (counter > 0)
        std::cerr << "Removed " << counter << " k-mers" << ::std::endl;

    return result;
}

template <>
struct FunctorComplement<AminoAcid>:
    public::std::function<AminoAcid(AminoAcid)>
{
    inline AminoAcid operator()(AminoAcid x) const
    {
        return x;
    }

};

} // namespace seqan

namespace stellar
{

namespace app
{

template <typename TSequence, typename TId>
void _mergeMatchesIntoFirst(StringSet<QueryMatches<StellarMatch<TSequence const, TId> > > & matches1,
                            StringSet<QueryMatches<StellarMatch<TSequence const, TId> > > & matches2)
{
    // merge matches and reverseMatches
    if (length(matches1) != length(matches2))
        throw std::runtime_error{"Matches mismatch"};

    for (size_t i = 0; i < length(matches1); ++i)
    {
        matches1[i].mergeIn(matches2[i]);
    }
}

template <typename TAlphabet, typename TId>
inline StellarOutputStatistics
_stellarOnWholeDatabase(StringSet<String<TAlphabet> > const & databases,
                        StringSet<TId> const & databaseIDs,
                        StringSet<String<TAlphabet> > const & queries,
                        StringSet<TId> const & queryIDs,
                        bool const databaseStrand,
                        StellarOptions & options,
                        StellarSwiftPattern<TAlphabet> & swiftPattern,
                        std::vector<size_t> & disabledQueryIDs,
                        std::ofstream & outputFile)
{
    using TSequence = String<TAlphabet>;

    // container for eps-matches
    StringSet<QueryMatches<StellarMatch<TSequence const, TId> > > matches;
    resize(matches, length(queries));

    StellarComputeStatisticsCollection computeStatistics{length(databases)};

    #pragma omp parallel default(none) firstprivate(databaseStrand) shared(databases, options, swiftPattern, databaseIDs, matches, computeStatistics)
    {
        StringSet<QueryMatches<StellarMatch<TSequence const, TId> > > localMatches;
        resize(localMatches, length(matches));

        StellarOptions localOptions = options;
        StellarComputeStatisticsPartialCollection localPartialStatistics{computeStatistics.size()};

        StellarSwiftPattern<TAlphabet> localSwiftPattern = swiftPattern;

        #pragma omp for nowait
        for (size_t i = 0; i < length(databases); ++i)
        {
            String<TAlphabet> const & database = databases[i];
            CharString const & databaseID = databaseIDs[i];
            size_t const databaseRecordID = i;

            auto getQueryMatches = [&](auto const & pattern) -> QueryMatches<StellarMatch<TSequence const, TId> > &
            {
                return value(localMatches, pattern.curSeqNo);
            };

            auto isPatternDisabled = [&](StellarSwiftPattern<TAlphabet> & pattern) -> bool {
                QueryMatches<StellarMatch<TSequence const, TId> > & queryMatches = getQueryMatches(pattern);
                return queryMatches.disabled;
            };

            auto onAlignmentResult = [&](auto & alignment) -> bool {
                QueryMatches<StellarMatch<TSequence const, TId> > & queryMatches = getQueryMatches(localSwiftPattern);

                StellarMatch<TSequence const, TId> match(alignment, databaseID, databaseStrand);
                length(match);  // DEBUG: Contains assertion on clipping.

                // success
                return _insertMatch(
                    queryMatches,
                    match,
                    localOptions.minLength,
                    localOptions.disableThresh,
                    // compactThresh is basically an output-parameter; will be updated in kernel and propagated back
                    // outside of this function, the reason why StellarOptions can't be passed as const to this function.
                    // TODO: We might want to make this tied to the QueryMatches itself, as it should know then to consolidate
                    // the matches
                    localOptions.compactThresh,
                    localOptions.numMatches
                );
            };

            StellarDatabaseSegment<TAlphabet> databaseSegment{database, 0u, length(database)};
            StellarComputeStatistics statistics
                = _stellarOnOne(databaseSegment, localSwiftPattern, localOptions, isPatternDisabled, onAlignmentResult);

            localPartialStatistics.updateByRecordID(databaseRecordID, statistics);
        }

        #pragma omp critical
        {
            // linear in database count
            computeStatistics.mergePartialIn(localPartialStatistics);
            // linear in queries count
            _mergeMatchesIntoFirst(matches, localMatches);
        }
    } // parallel region - end

    // standard output:
    std::cerr << std::endl; // swift filter output is on same line
    for (size_t i = 0; i < length(databases); ++i)
    {
        CharString const & databaseID = databaseIDs[i];
        std::cout << "  " << databaseID;
        if (!databaseStrand)
            std::cout << ", complement";
        std::cout << std::flush;

        if (options.verbose)
        {
            StellarComputeStatistics const & statistics = computeStatistics[i];
            _printStellarKernelStatistics(statistics);
        }
        std::cout << std::endl;
    }

    for (size_t queryID = 0; queryID < length(matches); ++queryID)
    {
        QueryMatches<StellarMatch<TSequence const, TId>> & queryMatches = value(matches, queryID);

        queryMatches.removeOverlapsAndCompactMatches(options.disableThresh,
                                                     /*compactThresh*/ 0,
                                                     options.minLength,
                                                     options.numMatches);

        if (queryMatches.disabled)
            disabledQueryIDs.push_back(queryID);
    }

    // if databaseStrand == true always outputs
    // if databaseStrand == false only outputs if TAlphabet == Dna5 or TAlphabet == Rna5
    if (databaseStrand || IsSameType<TAlphabet, Dna5>::VALUE || IsSameType<TAlphabet, Rna5>::VALUE)
    {
        // adjust length for each matches of a single query (only for dna5 and rna5)
        _postproccessLengthAdjustment(matches);

        // output matches on positive database strand
        _writeAllQueryMatchesToFile(matches, queryIDs, databaseStrand, options.outputFormat, outputFile);
    }

    return _computeOutputStatistics(matches);
}

///////////////////////////////////////////////////////////////////////////////
// Initializes a Pattern object with the query sequences,
//  and calls _stellarOnOne for each database sequence
template <typename TAlphabet, typename TId>
inline bool
_stellarOnAll(StringSet<String<TAlphabet>> & databases,
              StringSet<TId> const & databaseIDs,
              StringSet<String<TAlphabet>> const & queries,
              StringSet<TId> const & queryIDs,
              StellarOptions & options)
{
    // pattern
    StellarIndex<TAlphabet> stellarIndex{queries, options};
    StellarSwiftPattern<TAlphabet> swiftPattern = stellarIndex.createSwiftPattern();

    if (options.verbose)
        swiftPattern.params.printDots = true;

    // Construct index
    std::cout << "Constructing index..." << std::endl;
    stellarIndex.construct();
    std::cout << std::endl;

    std::cout << "Aligning all query sequences to database sequence..." << std::endl;

    std::ofstream outputFile(toCString(options.outputFile), ::std::ios_base::out | ::std::ios_base::app);
    if (!outputFile.is_open()) {
        std::cerr << "Could not open output file." << std::endl;
        return false;
    }

    std::vector<size_t> disabledQueryIDs{};

    // compute and print output statistics
    StellarOutputStatistics outputStatistics{};

    // positive database strand
    if (options.forward)
    {
        outputStatistics
            = _stellarOnWholeDatabase(databases,
                                      databaseIDs,
                                      queries,
                                      queryIDs,
                                      true,
                                      options,
                                      swiftPattern,
                                      disabledQueryIDs,
                                      outputFile);
    }

    // negative (reverse complemented) database strand
    bool const reverse = options.reverse && options.alphabet != "protein" && options.alphabet != "char";
    if (reverse)
    {
        for (size_t i = 0; i < length(databases); ++i)
            reverseComplement(databases[i]);

        StellarOutputStatistics statistics
            = _stellarOnWholeDatabase(databases,
                                      databaseIDs,
                                      queries,
                                      queryIDs,
                                      false,
                                      options,
                                      swiftPattern,
                                      disabledQueryIDs,
                                      outputFile);

        outputStatistics.mergeIn(statistics);
    }
    std::cout << std::endl;

    bool const writeDisabledQueriesFile = options.disableThresh != std::numeric_limits<unsigned>::max();

    // Writes disabled query sequences to disabledFile.
    if (writeDisabledQueriesFile)
    {
        std::ofstream disabledQueriesFile(toCString(options.disabledQueriesFile),
                                          ::std::ios_base::out | ::std::ios_base::app);

        if (!disabledQueriesFile.is_open()) {
            std::cerr << "Could not open file for disabled queries." << std::endl;
            return false;
        }

        // write disabled query file
        _writeDisabledQueriesToFastaFile(disabledQueryIDs, queryIDs, queries, disabledQueriesFile);
    }

    _writeOutputStatistics(outputStatistics, options.verbose, writeDisabledQueriesFile);

    return true;
}

template <typename TId>
inline bool
_checkUniqueId(std::set<TId> & uniqueIds, TId const & id)
{
    TId shortId;
    typedef typename Iterator<TId const>::Type TIterator;

    TIterator it = begin(id);
    TIterator itEnd = end(id);

    // (cut at first whitespace)
    while (it != itEnd && *it > ' ')
    {
        appendValue(shortId, *it);
        ++it;
    }

    if (uniqueIds.count(shortId) == 0)
    {
        uniqueIds.insert(shortId);
        return 1;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Imports sequences from a file,
//  stores them in the StringSet seqs and their identifiers in the StringSet ids
template <typename TSequence, typename TId>
inline bool
_importSequences(CharString const & fileName,
                 CharString const & name,
                 StringSet<TSequence> & seqs,
                 StringSet<TId> & ids)
{
    SeqFileIn inSeqs;
    if (!open(inSeqs, (toCString(fileName))))
    {
        std::cerr << "Failed to open " << name << " file." << std::endl;
        return false;
    }

    std::set<TId> uniqueIds; // set of short IDs (cut at first whitespace)
    bool idsUnique = true;

    TSequence seq;
    TId id;
    unsigned seqCount = 0;
    for (; !atEnd(inSeqs); ++seqCount)
    {
        readRecord(id, seq, inSeqs);

        idsUnique &= _checkUniqueId(uniqueIds, id);

        appendValue(seqs, seq, Generous());
        appendValue(ids, id, Generous());
    }

    std::cout << "Loaded " << seqCount << " " << name << " sequence" << ((seqCount > 1) ? "s." : ".") << std::endl;
    if (!idsUnique)
        std::cerr << "WARNING: Non-unique " << name << " ids. Output can be ambiguous.\n";
    return true;
}

///////////////////////////////////////////////////////////////////////////////
// Parses and outputs parameters, calls _stellarOnAll().
template <typename TAlphabet>
int mainWithOptions(StellarOptions & options, String<TAlphabet>)
{
    typedef String<TAlphabet> TSequence;

    // set threads
    omp_set_num_threads(options.threadCount);

    // output file names
    stellar::app::_writeFileNames(options);

    // output parameters
    stellar::app::_writeSpecifiedParams(options);
    stellar::app::_writeCalculatedParams(options);

    // import query sequences
    StringSet<TSequence> queries;
    StringSet<CharString> queryIDs;
    if (!_importSequences(options.queryFile, "query", queries, queryIDs))
        return 1;

    // import database sequence
    StringSet<TSequence> databases;
    StringSet<CharString> databaseIDs;
    if (!_importSequences(options.databaseFile, "database", databases, databaseIDs))
        return 1;

    std::cout << std::endl;
    stellar::app::_writeMoreCalculatedParams(options, databases, queries);

    // open output files
    std::ofstream file;
    file.open(toCString(options.outputFile));
    if (!file.is_open())
    {
        std::cerr << "Could not open output file." << std::endl;
        return 1;
    }
    file.close();

    if (options.disableThresh != std::numeric_limits<unsigned>::max())
    {
        std::ofstream daFile;
        daFile.open(toCString(options.disabledQueriesFile));
        if (!daFile.is_open())
        {
            std::cerr << "Could not open file for disabled queries." << std::endl;
            return 1;
        }
        daFile.close();
    }

    // stellar on all databases and queries writing results to file

    double startTime = sysTime();
    if (!_stellarOnAll(databases, databaseIDs, queries, queryIDs, options))
        return 1;

    if (options.verbose && options.noRT == false)
        std::cout << "Running time: " << sysTime() - startTime << "s" << std::endl;

    return 0;
}

} // namespace stellar::app

} // namespace stellar
