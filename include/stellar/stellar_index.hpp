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

#ifndef SEQAN_HEADER_STELLAR_INDEX_HPP
#define SEQAN_HEADER_STELLAR_INDEX_HPP

#include <seqan/index.h>

#include <seqan3/std/span>

#include <stellar/stellar_types.hpp> // StellarOptions

namespace stellar
{
using namespace seqan;

template <typename TAlphabet, typename TString = String<TAlphabet>, typename TInfixSegment = seqan::Segment<TString const, seqan::InfixSegment>>
using StellarQGramStringSet = StringSet<TInfixSegment, Owner<> >;

template <typename TAlphabet>
using StellarQGramIndex = Index<StellarQGramStringSet<TAlphabet> const, IndexQGram<SimpleShape, OpenAddressing> >;

template <typename TAlphabet>
using StellarSwiftPattern = Pattern<StellarQGramIndex<TAlphabet>, Swift<SwiftLocal> >;

template <typename TAlphabet>
using StellarSwiftFinder = Finder<Segment<String<TAlphabet> const, InfixSegment> const, Swift<SwiftLocal> >;

template <typename TAlphabet>
struct StellarIndex
{
    using TSequence = seqan::String<TAlphabet>;
    using TInfixSegment = seqan::Segment<seqan::String<TAlphabet> const, seqan::InfixSegment>;
    using TQGramStringSet = StellarQGramStringSet<TAlphabet>;

    template <typename TSpec>
    StellarIndex(StringSet<TSequence, TSpec> const & queries, StellarOptions const & options)
        : StellarIndex{convertImplStringSet(queries), options}
    {}

    StellarIndex(std::span<TInfixSegment> const & queries, StellarOptions const & options)
        : StellarIndex{convertImplSpan(queries), options}
    {}

    void construct()
    {
        indexRequire(qgramIndex, QGramSADir());
    }

    StellarSwiftPattern<TAlphabet> createSwiftPattern()
    {
        return {qgramIndex};
    }

    static StellarQGramIndex<TAlphabet> & qgramIndexFromPattern(StellarSwiftPattern<TAlphabet> & pattern)
    {
        return host(pattern);
    }

    static TQGramStringSet const & sequencesFromPattern(StellarSwiftPattern<TAlphabet> & pattern)
    {
        return sequencesFromQGramIndex(qgramIndexFromPattern(pattern));
    }

    static TQGramStringSet const & sequencesFromQGramIndex(StellarQGramIndex<TAlphabet> & index)
    {
        return indexText(index);
    }

    TQGramStringSet const dependentQueries;
    StellarQGramIndex<TAlphabet> qgramIndex;

private:
    template <typename TOtherQGramStringSet, typename = std::enable_if_t<std::is_same_v<TOtherQGramStringSet, TQGramStringSet>>>
    StellarIndex(TOtherQGramStringSet && queries, StellarOptions const & options)
        : dependentQueries{std::forward<TOtherQGramStringSet>(queries)}, qgramIndex{dependentQueries}
    {
        resize(indexShape(qgramIndex), options.qGram);
        cargo(qgramIndex).abundanceCut = options.qgramAbundanceCut;
    }

    template <typename TSpec>
    static StellarQGramStringSet<TAlphabet> convertImplStringSet(StringSet<TSequence, TSpec> const & queries)
    {
        StellarQGramStringSet<TAlphabet> dependentQueries;
        for (TSequence const & query: queries)
            seqan::appendValue(dependentQueries, seqan::infix(query, 0, seqan::length(query)));

        return dependentQueries;
    }

    static StellarQGramStringSet<TAlphabet> convertImplSpan(std::span<TInfixSegment> const & queries)
    {
        StellarQGramStringSet<TAlphabet> dependentQueries;
        for (TInfixSegment const & query: queries)
            seqan::appendValue(dependentQueries, query);

        return dependentQueries;
    }
};

} // namespace stellar

#endif
