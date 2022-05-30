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
#include <stellar/database_id_map.hpp>

#include <stellar/parallel/compute_statistics_collection.hpp>

#include <stellar/prefilter/no_query_prefilter.hpp>
#include <stellar/prefilter/whole_database_agent_splitter.hpp>

#include <stellar/app/stellar.diagnostics.hpp>

#include <stellar/app/stellar.diagnostics.tpp>

namespace stellar
{
using namespace seqan;
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
bool _shouldWriteOutputFile(bool const databaseStrand, StringSet<QueryMatches<StellarMatch<String<TAlphabet> const, TId> > > const & matches)
{
    // if databaseStrand == true always outputs
    // if databaseStrand == false only outputs if TAlphabet == Dna5 or TAlphabet == Rna5
    return databaseStrand || IsSameType<TAlphabet, Dna5>::VALUE || IsSameType<TAlphabet, Rna5>::VALUE;
}

template <typename TAlphabet, typename TId>
void _postproccessQueryMatches(bool const databaseStrand, StellarOptions const & options, StringSet<QueryMatches<StellarMatch<String<TAlphabet> const, TId> > > & matches, std::vector<size_t> & disabledQueryIDs)
{
    using TSequence = String<TAlphabet>;

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

    // adjust length for each matches of a single query (only for dna5 and rna5)
    // TODO: WHY? This seems like an arbitrary restriction
    if (_shouldWriteOutputFile(databaseStrand, matches))
        _postproccessLengthAdjustment(matches);
}

template <typename TAlphabet, typename TId>
inline StellarComputeStatisticsCollection
_parallelPrefilterStellar(
    StringSet<String<TAlphabet> > const & databases,
    StringSet<TId> const & databaseIDs,
    StringSet<String<TAlphabet> > const & queries,
    bool const databaseStrand,
    StellarOptions const & options,
    StellarSwiftPattern<TAlphabet> & swiftPattern,
    StringSet<QueryMatches<StellarMatch<String<TAlphabet> const, TId> > > & matches)
{
    using TSequence = String<TAlphabet>;

    using TQueryFilter = StellarSwiftPattern<TAlphabet>;
    using TSplitter = WholeDatabaseAgentSplitter;
    using TPrefilter = NoQueryPrefilter<TAlphabet, TQueryFilter, TSplitter>;
    using TPrefilterAgent = typename TPrefilter::Agent;
    using TDatabaseSegments = typename TPrefilter::TDatabaseSegments;

    static constexpr auto _verificationMethodVisit =
    [](StellarVerificationMethod verificationMethod, auto && visitor_fn)
        -> StellarComputeStatistics
    {
        if (verificationMethod == StellarVerificationMethod{AllLocal{}})
            return visitor_fn(AllLocal());
        else if (verificationMethod == StellarVerificationMethod{BestLocal{}})
            return visitor_fn(BestLocal());
        else if (verificationMethod == StellarVerificationMethod{BandedGlobal{}})
            return visitor_fn(BandedGlobal());
        else if (verificationMethod == StellarVerificationMethod{BandedGlobalExtend{}})
            return visitor_fn(BandedGlobalExtend());
        return StellarComputeStatistics{};
    };

    TPrefilter prefilter{databases, TQueryFilter{swiftPattern} /*copy pattern*/, TSplitter{}};

    StellarComputeStatisticsCollection computeStatistics{length(databases)};

    DatabaseIDMap<TAlphabet> databaseIDMap{databases, databaseIDs};

    std::vector<TPrefilterAgent> prefilterAgents = prefilter.agents(options.threadCount, options.minLength);

    #pragma omp parallel for num_threads(prefilterAgents.size()) default(none) firstprivate(databaseStrand) shared(std::cout, prefilterAgents, options, matches, databaseIDMap, computeStatistics)
    for (TPrefilterAgent & agent: prefilterAgents)
    {
        StringSet<QueryMatches<StellarMatch<TSequence const, TId> > > localMatches;
        resize(localMatches, length(matches));

        StellarOptions localOptions = options;
        StellarComputeStatisticsPartialCollection localPartialStatistics{computeStatistics.size()};

        agent.prefilter([&](TDatabaseSegments const & databaseSegments, TQueryFilter localSwiftPattern)
        {
            for (StellarDatabaseSegment<TAlphabet> const & databaseSegment : databaseSegments)
            {
                String<TAlphabet> const & database = databaseSegment.underlyingDatabase();
                size_t const databaseRecordID = databaseIDMap.recordID(database);

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

                    TId const & databaseID = databaseIDMap.databaseIDs[databaseRecordID];
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

                // finder
                StellarSwiftFinder<TAlphabet> swiftFinder(databaseSegment.asInfixSegment(), localOptions.minRepeatLength, localOptions.maxRepeatPeriod);

                StellarComputeStatistics statistics = _verificationMethodVisit(
                    localOptions.verificationMethod,
                    [&](auto tag) -> StellarComputeStatistics
                    {
                        using TTag = decltype(tag);
                        SwiftHitVerifier<TTag> swiftVerifier
                        {
                            STELLAR_DESIGNATED_INITIALIZER(.epsilon = , localOptions.epsilon),
                            STELLAR_DESIGNATED_INITIALIZER(.minLength = , localOptions.minLength),
                            STELLAR_DESIGNATED_INITIALIZER(.xDrop = , localOptions.xDrop)
                        };

                        return _stellarKernel(swiftFinder, localSwiftPattern, swiftVerifier, isPatternDisabled, onAlignmentResult);
                    });

                localPartialStatistics.updateByRecordID(databaseRecordID, statistics);
            }
        });

        #pragma omp critical
        {
            // linear in database count
            computeStatistics.mergePartialIn(localPartialStatistics);
            // linear in queries count
            _mergeMatchesIntoFirst(matches, localMatches);
        }
    } // parallel region - end

    return computeStatistics;
}

///////////////////////////////////////////////////////////////////////////////
// Initializes a Pattern object with the query sequences,
//  and calls _parallelPrefilterStellar for each database sequence
template <typename TAlphabet, typename TId>
inline bool
_stellarMain(
    StringSet<String<TAlphabet>> & databases,
    StringSet<TId> const & databaseIDs,
    StringSet<String<TAlphabet>> const & queries,
    StringSet<TId> const & queryIDs,
    StellarOptions const & options,
    std::ofstream & outputFile,
    std::ofstream & disabledQueriesFile)
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

    std::vector<size_t> disabledQueryIDs{};

    // compute and print output statistics
    StellarOutputStatistics outputStatistics{};

    // positive database strand
    if (options.forward)
    {
        // container for eps-matches
        StringSet<QueryMatches<StellarMatch<String<TAlphabet> const, TId> > > forwardMatches;
        resize(forwardMatches, length(queries));

        constexpr bool databaseStrand = true;

        StellarComputeStatisticsCollection computeStatistics =
            _parallelPrefilterStellar(
                databases,
                databaseIDs,
                queries,
                databaseStrand,
                options,
                swiftPattern,
                forwardMatches);

        // standard output:
        _printParallelPrefilterStellarStatistics(options.verbose, databaseStrand, databaseIDs, computeStatistics);

        _postproccessQueryMatches(databaseStrand, options, forwardMatches, disabledQueryIDs);

        if (_shouldWriteOutputFile(databaseStrand, forwardMatches))
        {
            // output forwardMatches on positive database strand
            _writeAllQueryMatchesToFile(forwardMatches, queryIDs, databaseStrand, options.outputFormat, outputFile);
        }

        outputStatistics = _computeOutputStatistics(forwardMatches);
    }

    // negative (reverse complemented) database strand
    bool const reverse = options.reverse && options.alphabet != "protein" && options.alphabet != "char";
    if (reverse)
    {
        for (size_t i = 0; i < length(databases); ++i)
            reverseComplement(databases[i]);

        // container for eps-matches
        StringSet<QueryMatches<StellarMatch<String<TAlphabet> const, TId> > > reverseMatches;
        resize(reverseMatches, length(queries));

        constexpr bool databaseStrand = false;

        StellarComputeStatisticsCollection computeStatistics =
            _parallelPrefilterStellar(
                databases,
                databaseIDs,
                queries,
                databaseStrand,
                options,
                swiftPattern,
                reverseMatches);

        // standard output:
        _printParallelPrefilterStellarStatistics(options.verbose, databaseStrand, databaseIDs, computeStatistics);

        _postproccessQueryMatches(databaseStrand, options, reverseMatches, disabledQueryIDs);

        if (_shouldWriteOutputFile(databaseStrand, reverseMatches))
        {
            // output reverseMatches on negative database strand
            _writeAllQueryMatchesToFile(reverseMatches, queryIDs, databaseStrand, options.outputFormat, outputFile);
        }

        outputStatistics.mergeIn(_computeOutputStatistics(reverseMatches));
    }
    std::cout << std::endl;

    // Writes disabled query sequences to disabledFile.
    if (disabledQueriesFile.is_open())
    {
        // write disabled query file
        _writeDisabledQueriesToFastaFile(disabledQueryIDs, queryIDs, queries, disabledQueriesFile);
    }

    _writeOutputStatistics(outputStatistics, options.verbose, disabledQueriesFile.is_open());

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
// Parses and outputs parameters, calls _stellarMain().
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
    std::ofstream outputFile(toCString(options.outputFile), ::std::ios_base::out | ::std::ios_base::app);
    if (!outputFile.is_open())
    {
        std::cerr << "Could not open output file." << std::endl;
        return 1;
    }

    std::ofstream disabledQueriesFile;
    if (options.disableThresh != std::numeric_limits<unsigned>::max())
    {
        disabledQueriesFile.open(toCString(options.disabledQueriesFile), ::std::ios_base::out | ::std::ios_base::app);
        if (!disabledQueriesFile.is_open())
        {
            std::cerr << "Could not open file for disabled queries." << std::endl;
            return 1;
        }
    }

    // stellar on all databases and queries writing results to file

    double startTime = sysTime();
    if (!_stellarMain(databases, databaseIDs, queries, queryIDs, options, outputFile, disabledQueriesFile))
        return 1;

    if (options.verbose && options.noRT == false)
        std::cout << "Running time: " << sysTime() - startTime << "s" << std::endl;

    return 0;
}

} // namespace stellar::app

} // namespace stellar
