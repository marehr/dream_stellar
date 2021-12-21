
#pragma once

#include "../stellar_types.h"

namespace stellar
{

namespace app
{

///////////////////////////////////////////////////////////////////////////////
// Calculates parameters from parameters in options object and writes them to std::cout
void _writeCalculatedParams(StellarOptions & options);

///////////////////////////////////////////////////////////////////////////////
// Writes user specified parameters from options object to std::cout
void _writeSpecifiedParams(StellarOptions const & options);

///////////////////////////////////////////////////////////////////////////////
// Writes file name from options object to std::cout
void _writeFileNames(StellarOptions const & options);

///////////////////////////////////////////////////////////////////////////////
// Calculates parameters from parameters in options object and from sequences and writes them to std::cout
template <typename TStringSet>
void _writeMoreCalculatedParams(StellarOptions const & options, TStringSet const & databases, TStringSet const & queries);

void _writeOutputStatistics(StellarOutputStatistics const & statistics, bool const verbose, bool const writeDisabledQueriesFile);

void _printStellarKernelStatistics(StellarComputeStatistics const & statistics);

} // namespace stellar::app

} // namespace stellar
