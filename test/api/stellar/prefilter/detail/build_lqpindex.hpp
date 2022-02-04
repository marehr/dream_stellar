#pragma once

#include <stellar/prefilter/detail/valik_local_query_prefilter_index.hpp>

auto build_lqpindex()
{
    stellar::detail::valik_local_query_prefilter_index lqpindex
    {
        // ibf parameters
        ._ibf {
            seqan3::bin_count{5u /*arguments->bins*/},
            seqan3::bin_size{4096u / 1u/*arguments->bits / arguments->parts*/},
            seqan3::hash_function_count{2u/*arguments->hash*/}
        },

        // hash parameters
        ._kmer_size{seqan3::ungapped{4u/*arguments->kmer_size*/}},
        ._window_size{seqan3::window_size{6u/*arguments->window_size*/}}
    };

    return lqpindex;
}

auto build_lqpindex2()
{
    stellar::detail::valik_local_query_prefilter_index lqpindex
    {
        // ibf parameters
        ._ibf {
            seqan3::bin_count{5u},
            seqan3::bin_size{524'288'000u / 16u}, /* = 32'768'000u*/
            // 32MB for 1MB sequence file LOL
            seqan3::hash_function_count{2u/*arguments->hash*/}
        },

        // hash parameters
        ._kmer_size{seqan3::ungapped{19u/*arguments->kmer_size*/}},
        ._window_size{seqan3::window_size{24u/*arguments->window_size*/}}
    };

    return lqpindex;
}
