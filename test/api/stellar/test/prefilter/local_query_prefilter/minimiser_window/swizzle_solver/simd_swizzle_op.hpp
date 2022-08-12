
#pragma once

#include <cassert>
#include <variant>

#include "simd_int32x.hpp"

struct argument_m256i
{
    int position{-1};
};

struct argument_constant_m256i
{
    simd_int32x8 constant;
};

struct argument_any_m256i
{
    simd_int32x8 constant;
};

struct argument_imm8
{
    uint8_t constant;
};

struct argument_none
{};

template <typename argument0_t, typename argument1_t, typename argument2_t>
struct simd_swizzle_op
{
    simd_int32x8 shuffle{};

    std::string_view name{};

    argument0_t arg0{};
    argument1_t arg1{};
    argument2_t arg2{};

    argument_m256i return_argument{};

    static constexpr std::size_t arg_count =
        std::is_same_v<argument0_t, argument_m256i> +
        std::is_same_v<argument1_t, argument_m256i> +
        std::is_same_v<argument2_t, argument_m256i>;

    simd_int32x8 operator()(simd_int32x8 const & memory) const requires (arg_count == 1)
    {
        simd_int32x8 data{};

        for (int i = 0; i < 8; ++i)
        {
            assert(shuffle[i] < 8);
            data[i] = memory[shuffle[i]];
        }

        return data;
    }

    simd_int32x8 operator()(simd_int32x8 const & memory1, simd_int32x8 const & memory2) const requires (arg_count == 2)
    {
        simd_int32x8 data{};

        for (int i = 0; i < 8; ++i)
        {
            assert(shuffle[i] < 16);
            if (shuffle[i] < 8)
                data[i] = memory1[shuffle[i]];
            else
                data[i] = memory2[shuffle[i] - 8];
        }

        return data;
    }

    template <typename char_t>
    friend std::basic_ostream<char_t> & operator<<(std::basic_ostream<char_t> & cout, simd_swizzle_op const & swizzle_op)
    {
        auto print_argument = [&](auto && argument)
        {
            using argument_t = std::decay_t<decltype(argument)>;
            if constexpr (std::is_same_v<argument_t, argument_m256i>)
            {
                cout << "_m256i[" << (argument.position == -1 ? "?" : std::to_string(argument.position)) << "]";
            }
            else if constexpr (std::is_same_v<argument_t, argument_constant_m256i>)
            {
                cout << "_m256i{" << argument.constant << "}";
            }
            else if constexpr (std::is_same_v<argument_t, argument_any_m256i>)
            {
                cout << "_m256i{ANY}";
            }
            else if constexpr (std::is_same_v<argument_t, argument_imm8>)
            {
                cout << "imm8{" << (int)argument.constant << "}";
            }
            else if constexpr (std::is_same_v<argument_t, argument_none>)
            {
                cout << "Q";
            }
        };

        print_argument(swizzle_op.return_argument);
        cout << " = ";
        cout << swizzle_op.name;
        cout << "(";
        print_argument(swizzle_op.arg0);
        cout << ", ";
        print_argument(swizzle_op.arg1);
        cout << ", ";
        print_argument(swizzle_op.arg2);
        cout << ")";
        return cout;
    }
};

using simd_swizzle_opv = std::variant<
    simd_swizzle_op<argument_m256i, argument_any_m256i, argument_imm8>,
    simd_swizzle_op<argument_m256i, argument_m256i, argument_imm8>,
    simd_swizzle_op<argument_m256i, argument_imm8, argument_none>,
    simd_swizzle_op<argument_m256i, argument_m256i, argument_none>//,
>;
