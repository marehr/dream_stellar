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

#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

#include "stellar.h"
#include "stellar_output.h"

#include "stellar/stellar.arg_parser.hpp"
#include "stellar/stellar.diagnostics.hpp"
#include "stellar/stellar.main.hpp"

#ifndef STELLAR_PARALLEL_BUILD
#include "stellar/stellar.arg_parser.cpp"
#include "stellar/stellar.diagnostics.cpp"
#include "stellar/stellar.main.tpp"
#endif // STELLAR_PARALLEL_BUILD

// TODO(holtgrew): Move this into a SeqAn misc module.

// not supported anymore in vc2015
// https://msdn.microsoft.com/en-us/library/bb531344.aspx
class ScientificNotationExponentOutputNormalizer
{
public:
    unsigned _oldExponentFormat;

    ScientificNotationExponentOutputNormalizer() :
        _oldExponentFormat(0)
    {
    }

    ~ScientificNotationExponentOutputNormalizer()
    {
    }
};

///////////////////////////////////////////////////////////////////////////////
// Program entry point.
int main(int argc, const char * argv[])
{
    // Makes sure that printing doubles in scientific notation is normalized on all platforms.
    ScientificNotationExponentOutputNormalizer scientificNotationNormalizer;

    // command line parsing
    seqan::ArgumentParser parser("stellar");

    stellar::StellarOptions options{};
    stellar::app::_setParser(parser);
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    if (res == seqan::ArgumentParser::PARSE_OK)
        res = stellar::app::_parseOptions(parser, options);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    if (options.alphabet == "dna")
        stellar::app::mainWithOptions(options, seqan::String<seqan::Dna>());
    else if (options.alphabet == "dna5")
        stellar::app::mainWithOptions(options, seqan::String<seqan::Dna5>());
    else if (options.alphabet == "rna")
        stellar::app::mainWithOptions(options, seqan::String<seqan::Rna>());
    else if (options.alphabet == "rna5")
        stellar::app::mainWithOptions(options, seqan::String<seqan::Rna5>());
    else if (options.alphabet == "protein")
        stellar::app::mainWithOptions(options, seqan::String<seqan::AminoAcid>());
    else if (options.alphabet == "char")
        stellar::app::mainWithOptions(options, seqan::String<char>());
}
