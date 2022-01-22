
#pragma once

#include <functional>

#include <seqan/basic.h>
#include <seqan/modifier/modifier_functors.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>

// class SimpleType
// template <>
template <typename TValue, typename TSpec>
struct seqan3::custom::alphabet<seqan::SimpleType<TValue, TSpec>>
{
    static_assert(std::is_same_v<TValue, uint8_t>);

    using alphabet_t = seqan::SimpleType<TValue, TSpec>;
    using rank_t = uint8_t;
    using char_t = char;

    static constexpr size_t alphabet_size = seqan::ValueSize<alphabet_t>::VALUE;

    static /*constexpr*/ rank_t to_rank(alphabet_t const a) noexcept
    {
        return seqan::ordValue(a);
    }

    static /*constexpr*/ alphabet_t & assign_rank_to(rank_t const r, alphabet_t & a) noexcept
    {
        a = r;
        return a;
    }

    static /*constexpr*/ char_t to_char(alphabet_t const a) noexcept
    {
        return (char_t) a;
    }

    static /*constexpr*/ alphabet_t & assign_char_to(char_t const c, alphabet_t & a) noexcept
    {
        a = c;
        return a;
    }

    static /*constexpr*/ alphabet_t complement(alphabet_t const a) noexcept
    {
        return seqan::FunctorComplement<alphabet_t>{}(a);
    }
};
