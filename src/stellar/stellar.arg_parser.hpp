
#pragma once

#include <seqan/arg_parse.h>

#include "../stellar_types.h"

namespace stellar
{

namespace app
{

ArgumentParser::ParseResult
_parseOptions(ArgumentParser const & parser, StellarOptions & options);

void _setParser(ArgumentParser & parser);

} // namespace stellar::app

} // namespace stellar
