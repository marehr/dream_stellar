#include "stellar.diagnostics.hpp"

namespace stellar
{

namespace app
{

///////////////////////////////////////////////////////////////////////////////
// Calculates parameters from parameters in options object and writes them to std::cout
void _writeCalculatedParams(StellarOptions & options)
{
//IOREV _notio_
    StellarStatistics statistics{options};

    std::cout << "Calculated parameters:" << std::endl;
    if (statistics.kMerComputed)
    {
        options.qGram = (unsigned)statistics.kMerLength;
        std::cout << "  k-mer length : " << statistics.kMerLength << std::endl;
    }

    std::cout << "  s^min        : " << statistics.smin << std::endl;
    std::cout << "  threshold    : " << statistics.threshold << std::endl;
    std::cout << "  distance cut : " << statistics.distanceCut << std::endl;
    std::cout << "  delta        : " << statistics.delta << std::endl;
    std::cout << "  overlap      : " << statistics.overlap << std::endl;
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Writes user specified parameters from options object to std::cout
void _writeSpecifiedParams(StellarOptions const & options)
{
//IOREV _notio_
    // Output user specified parameters
    std::cout << "User specified parameters:" << std::endl;
    std::cout << "  minimal match length             : " << options.minLength << std::endl;
    std::cout << "  maximal error rate (epsilon)     : " << options.epsilon << std::endl;
    std::cout << "  maximal x-drop                   : " << options.xDrop << std::endl;
    if (options.qGram != (unsigned)-1)
        std::cout << "  k-mer (q-gram) length            : " << options.qGram << std::endl;
    std::cout << "  search forward strand            : " << ((options.forward) ? "yes" : "no") << std::endl;
    std::cout << "  search reverse complement        : " << ((options.reverse) ? "yes" : "no") << std::endl;
    std::cout << std::endl;

    std::cout << "  verification strategy            : " << options.fastOption << std::endl;
    if (options.disableThresh != (unsigned)-1)
    {
        std::cout << "  disable queries with more than   : " << options.disableThresh << " matches" << std::endl;
    }
    std::cout << "  maximal number of matches        : " << options.numMatches << std::endl;
    std::cout << "  duplicate removal every          : " << options.compactThresh << std::endl;
    if (options.maxRepeatPeriod != 1 || options.minRepeatLength != 1000)
    {
        std::cout << "  max low complexity repeat period : " << options.maxRepeatPeriod << std::endl;
        std::cout << "  min low complexity repeat length : " << options.minRepeatLength << std::endl;
    }
    if (options.qgramAbundanceCut != 1)
    {
        std::cout << "  q-gram abundance cut ratio       : " << options.qgramAbundanceCut << std::endl;
    }
    std::cout << "  threads                          : " << options.threadCount << std::endl;
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Writes file name from options object to std::cout
void _writeFileNames(StellarOptions const & options)
{
//IOREV _notio_
    std::cout << "I/O options:" << std::endl;
    std::cout << "  database file   : " << options.databaseFile << std::endl;
    std::cout << "  query file      : " << options.queryFile << std::endl;
    std::cout << "  alphabet        : " << options.alphabet << std::endl;
    std::cout << "  output file     : " << options.outputFile << std::endl;
    std::cout << "  output format   : " << options.outputFormat << std::endl;
    if (options.disableThresh != (unsigned)-1)
    {
        std::cout << "  disabled queries: " << options.disabledQueriesFile << std::endl;
    }
    std::cout << std::endl;
}

void _writeOutputStatistics(StellarOutputStatistics const & statistics, bool const verbose, bool const writeDisabledQueriesFile)
{
    std::cout << "# Eps-matches     : " << statistics.numMatches << std::endl;
    if (verbose) {
        if (statistics.numMatches > 0) {
            std::cout << "Longest eps-match : " << statistics.maxLength << std::endl;
            std::cout << "Avg match length  : " << statistics.totalLength / statistics.numMatches << std::endl;
        }
        if (writeDisabledQueriesFile)
            std::cout << "# Disabled queries: " << statistics.numDisabled << std::endl;
    }
}

void _printStellarKernelStatistics(StellarComputeStatistics const & statistics)
{
    if (statistics.numSwiftHits == 0)
        return;

    std::cout << std::endl << "    # SWIFT hits      : " << statistics.numSwiftHits;
    std::cout << std::endl << "    Longest hit       : " << statistics.maxLength;
    std::cout << std::endl << "    Avg hit length    : " << statistics.totalLength/statistics.numSwiftHits;
}

} // namespace stellar::app

} // namespace stellar
