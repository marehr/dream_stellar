#include <gtest/gtest.h>

// DEBUG: remove me
#include <seqan3/core/debug_stream.hpp>

#include "detail/build_lqpindex.hpp"

#include <stellar/stellar_types.hpp>
#include <stellar/prefilter/local_query_prefilter.hpp>

TEST(LocalQueryPrefilter, LocalQueryPrefilter)
{
    using TAlphabet = seqan::Dna5;

    seqan::StringSet<seqan::String<TAlphabet> > databases{};
    appendValue(databases, "TAGAC" "TAACCGCAGAACACG" "ACTCCTCTAC");
    appendValue(databases, "TCAC" "GCCAGTTTAATTTTA" "GGGCTCCAAGT");
    appendValue(databases, "TACGCATATCTGG" "TAACCGCAGAACACG" "AA");

    auto lqpindex = build_lqpindex();
    for (auto const & database : databases)
    {
        size_t const i = &database - &databases[0];
        lqpindex.insert_sequence(seqan3::bin_index{i}, database);
    }

    seqan::StringSet<seqan::String<TAlphabet> > queries{};
    appendValue(queries, "TAGCCG" "TAACCGCAGAACACG" "GACC"); // matches database[0,2]
    appendValue(queries, "TA" "GCCAGTTTAATTTTA" "TACCAAGA"); // matches database[1]
    appendValue(queries, "ACAGCTACCGT"); // matches nothing (does not match minLength)
    appendValue(queries, "TTGAGACGAAT" "AGAACACGACTCCTCTAC" "CTTGTGGCA"); // matches database[0]

    // Note: TQueryFilter are the queries themselve to check whith which queries the actual filter would have been
    //       constructed.
    using TPrefilter = stellar::LocalQueryPrefilter<TAlphabet>;
    using TPrefilterAgent = stellar::prefilter_agent<TAlphabet>;

    stellar::StellarOptions options{};
    options.minLength = 15u;
    options.threadCount = 1u;
    options.epsilon = {0, 1}; // 0.0

    TPrefilter prefilter{options, databases, queries, std::move(lqpindex)};

    using TDatabaseSegment = stellar::StellarDatabaseSegment<TAlphabet>;
    using TQuerySegment = seqan::Segment<seqan::String<TAlphabet> const, seqan::InfixSegment>;

    std::vector<size_t> seenQueryIDs{};

    for (TPrefilterAgent & prefilterAgent: prefilter.agents(options.threadCount, options))
    {
        using TDatabaseSegments = typename TPrefilter::TDatabaseSegments;
        using TQueryFilter = typename TPrefilter::TQueryFilter;

        prefilterAgent.prefilter([&](TDatabaseSegments databaseSegmentsSpan, TQueryFilter & queryFilter)
        {
            auto & sequences = stellar::StellarIndex<TAlphabet>::sequencesFromPattern(queryFilter);
            EXPECT_EQ(length(sequences), 1u);

            TQuerySegment const & querySegment = sequences[0];

            auto getQueryID = [&](TQuerySegment const & querySegment)
            {
                // NOTE: TQuerySegment, host of segment is the underlying query sequences.
                return &host(querySegment) - &queries[0];
            };
            auto getDatabaseID = [&](TDatabaseSegment const & databaseSegment)
            {
                return &databaseSegment.underlyingDatabase() - &databases[0];
            };
            auto getDatabaseIDs = [&](TDatabaseSegments const & databaseSegments)
            {
                std::vector<size_t> databaseIDs{};

                for (auto && databaseSegment : databaseSegments)
                    databaseIDs.push_back(getDatabaseID(databaseSegment));

                return databaseIDs;
            };

            std::vector<TDatabaseSegment> databaseSegments(databaseSegmentsSpan.begin(), databaseSegmentsSpan.end());
            std::sort(databaseSegments.begin(), databaseSegments.end());

            size_t const queryID = getQueryID(querySegment);

            EXPECT_LT(queryID, seqan::length(queries));

            seenQueryIDs.push_back(queryID);

            switch (queryID)
            {
                case 0:
                    EXPECT_EQ(querySegment, TQuerySegment(queries[0], 0, seqan::length(queries[0])));

                    EXPECT_EQ(databaseSegments.size(), 2);
                    EXPECT_EQ(getDatabaseIDs(databaseSegments), (std::vector<size_t>{0, 2}));
                    EXPECT_EQ(databaseSegments,
                        (std::vector<TDatabaseSegment>{
                            // D0: TAGAC [TAACCGCAGAACACG]ACTCCTCTAC
                            //            |||||||||||||||
                            // Q0: TAGCCG[TAACCGCAGAACACG]GACC"
                            {databases[0], 0, seqan::length(databases[0])},
                            // D2: TACGCATATCTG[GTAACCGCAGAACACG]AA
                            //                  ||||||||||||||||
                            // Q0: TAGCC       [GTAACCGCAGAACACG]GACC
                            {databases[2], 0, seqan::length(databases[2])}
                        })
                    );
                break;
                case 1:
                    EXPECT_EQ(querySegment, TQuerySegment(queries[1], 0, seqan::length(queries[1])));

                    EXPECT_EQ(databaseSegments.size(), 1);
                    EXPECT_EQ(getDatabaseIDs(databaseSegments), (std::vector<size_t>{1}));
                    EXPECT_EQ(databaseSegments,
                        (std::vector<TDatabaseSegment>{
                            // D1: TCAC[GCCAGTTTAATTTTA]GGGCTCCAAGT
                            //          |||||||||||||||
                            // Q1: TA  [GCCAGTTTAATTTTA]TACCAAGA
                            {databases[1], 0, seqan::length(databases[1])}
                        })
                    );
                break;
                // case 2: will been filtered out due to being to short
                case 3:
                    EXPECT_EQ(querySegment, TQuerySegment(queries[3], 0, seqan::length(queries[3])));

                    EXPECT_EQ(databaseSegments.size(), 1);
                    EXPECT_EQ(getDatabaseIDs(databaseSegments), (std::vector<size_t>{0}));
                    EXPECT_EQ(databaseSegments,
                        (std::vector<TDatabaseSegment>{
                            // D0: TAGACTAACCGC[AGAACACGACTCCTCTAC]
                            //                  ||||||||||||||||||
                            // Q3: TTGAGACGAAT [AGAACACGACTCCTCTAC]CTTGTGGCA
                            {databases[0], 0, seqan::length(databases[0])},
                            // // NOTE: This is a false-positive
                            // // D2: TACGCATATCTGGTAACCGC[AGAACACGA]A
                            // //                          |||||||||
                            // // Q3: TTGAGACGAAT         [AGAACACGA]CTCCTCTACCTTGTGGCA
                            // {databases[2], 0, seqan::length(databases[2])},
                        })
                    );
                break;
                default:
                    ADD_FAILURE() << "There should be only 4 queryIDs";
                break;
            };

#if STELLAR_DEBUG_MINIMISER
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << "queryID: " << queryID << std::endl;
            std::cout << "querySegment: " << queryFilter << std::endl;

            // collect all seen databases
            for (TDatabaseSegment databaseSegment : databaseSegments)
            {
                std::cout << "\tdatabaseID: " << getDatabaseID(databaseSegment) << std::endl;
                std::cout << "\tdatabaseSegment: " << databaseSegment.asInfixSegment() << std::endl;
            }

            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            std::cout << std::endl;
#endif // STELLAR_DEBUG_MINIMISER
        });
    }

    // Note: every query should have been seen completely
    std::sort(seenQueryIDs.begin(), seenQueryIDs.end());
    EXPECT_EQ(seenQueryIDs, (std::vector<size_t>
    {
        0, 1/*, 2 (will been filtered out due to beeing to short)*/, 3
    }));
}
