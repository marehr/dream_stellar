
#pragma once

#include "../stellar_types.h"

namespace stellar
{

namespace app
{

///////////////////////////////////////////////////////////////////////////////
// Parses and outputs parameters, calls _stellarOnAll().
template <typename TAlphabet>
int mainWithOptions(StellarOptions & options, String<TAlphabet>);

} // namespace stellar::app

} // namespace stellar
