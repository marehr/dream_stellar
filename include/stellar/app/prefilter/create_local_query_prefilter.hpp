
#pragma once

#include <stellar/stellar_types.hpp>

#include <stellar/app/prefilter/create_prefilter.hpp>

#include <stellar/prefilter/local_query_prefilter.hpp>

namespace stellar
{

namespace app
{

template <typename TAlphabet>
struct CreatePrefilter<LocalQueryPrefilter<TAlphabet>>
{
    static std::unique_ptr<stellar::prefilter<TAlphabet>>
    create(
        StellarOptions const & options,
        StringSet<String<TAlphabet> > const & databases,
        StringSet<String<TAlphabet> > const & queries,
        StellarSwiftPattern<TAlphabet> & swiftPattern)
    {
        using TPrefilter = LocalQueryPrefilter<TAlphabet>;

#if 0
        read length: 1000

        read: 100 -> bin 1
        read: 800 -> bin 2

        arguments.kmer_size: 19
        arguments.window_size: 24
        arguments.pattern_size: 50
        arguments.overlap: 49

    | bin_size      | minimiser | matched minimiser | matching pattern |
    | ------------- | --------- | ----------------- | ---------------- |
    | 524'288u      | 277599    | 119970            | 34552            |
    | 1'024'000u    | 277599    | 50898             | 2954             |
    | 2'048'000u    | 277599    | 18263             | 1938             |
    | 4'096'000u    | 277599    | 7084              | 1788             |
    | 8'192'000u    | 277599    | 3804              | 1747             |
    | 16'384'000u   | 277599    | 2933              | 1742             |
    | 32'768'000u   | 277599    | 2727              | 1741             |
    | 65'536'000u   | 277599    | 2660              | 1741             |
    | 131'072'000u  | 277599    | 2642              | 1741             |
    | 262'144'000u  | 277599    | 2639              | 1741             |
    | 524'288'000u  | 277599    | 2639              | 1741             |
    | 2'097'152'000 | 277599    | 2638              | 1741             |

#endif

        auto ibf_start = std::chrono::high_resolution_clock::now();
#if STELLAR_DEBUG_LOCAL_PREFILTER
        std::cout << "start IBF construction" << std::endl;
#endif // STELLAR_DEBUG_LOCAL_PREFILTER
        stellar::detail::valik_local_query_prefilter_index queryPrefilterIndex{
            // ibf parameters
            ._ibf {
                seqan3::bin_count{length(databases)},
                seqan3::bin_size{524'288'000u / 16u}, /* = 32'768'000u*/
                // 32MB for 1MB sequence file LOL
                seqan3::hash_function_count{2u/*arguments->hash*/}
            },

            // hash parameters
            ._kmer_size{seqan3::ungapped{19u/*arguments->kmer_size*/}},
            ._window_size{seqan3::window_size{24u/*arguments->window_size*/}}
        };

        if constexpr (!std::is_same_v<TAlphabet, char>)
        {
            for (size_t bin_id = 0; bin_id < length(databases); ++bin_id)
            {
                auto hash_view = queryPrefilterIndex.hash_view(databases[bin_id]);
#if STELLAR_DEBUG_LOCAL_PREFILTER
                std::cout << "databases[" << bin_id << "] has " << std::ranges::distance(hash_view) << " minimiser" << std::endl;
#endif // STELLAR_DEBUG_LOCAL_PREFILTER

                queryPrefilterIndex.insert_sequence(seqan3::bin_index{bin_id}, databases[bin_id]);
            }
        }

        TPrefilter prefilter{options, databases, queries, std::move(queryPrefilterIndex)};
        auto ibf_end = std::chrono::high_resolution_clock::now();
        double ibf_time = std::chrono::duration_cast<std::chrono::duration<double>>(ibf_end - ibf_start).count();

#if STELLAR_DEBUG_LOCAL_PREFILTER
        std::cout << "construct IBF in " << ibf_time << "s" << std::endl;
#endif // STELLAR_DEBUG_LOCAL_PREFILTER

        return std::make_unique<TPrefilter>(std::move(prefilter));
    }
};

} // namespace app

} // namespace stellar
