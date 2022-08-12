#include <algorithm>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <string_view>
#include <span>
#include <variant>
#include <vector>

void print_simd(std::basic_ostream<char> & cout, std::array<int, 8> const & simd)
{
    cout << "[";
    for (int i = 0; i < 8; ++i)
    {
        if (i != 0) cout << ",";
        cout << simd[i];
    }
    cout << "]";
}

struct simd_int32x8 : public std::array<int, 8>
{
    template <typename char_t>
    friend std::basic_ostream<char_t> & operator<<(std::basic_ostream<char_t> & cout, simd_int32x8 const & simd)
    {
        print_simd(cout, simd);
        return cout;
    }
};

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
    simd_swizzle_op<argument_m256i, argument_m256i, argument_imm8>//,
>;

struct simd_swizzle_state_rules
{
    simd_swizzle_state_rules() = default;

    simd_swizzle_state_rules(simd_swizzle_state_rules const &) = default;

    template <typename ...args_t>
        requires (sizeof...(args_t) > 1)
    simd_swizzle_state_rules(args_t && ...args)
        : transitions{simd_swizzle_opv{std::forward<args_t>(args)}...}
    {}

    std::vector<simd_swizzle_opv> transitions;
};

struct simd_swizzle_op_trace
{
    size_t size() const
    {
        return ops.size();
    }

    std::vector<simd_swizzle_opv> ops;

    template <typename char_t>
    friend std::basic_ostream<char_t> & operator<<(std::basic_ostream<char_t> & cout, simd_swizzle_op_trace const & trace)
    {
        cout << "[";
        size_t i = 0;
        for (simd_swizzle_opv const & opv : trace.ops)
        {
            std::visit([&cout](auto const & op)
            {
                cout << op;
            }, opv);

            if (++i < trace.ops.size())
                cout << " -> ";
        }
        cout << "]";
        return cout;
    }
};

template <size_t count>
struct simd_vars : public std::array<simd_int32x8, count>
{
    static constexpr simd_vars iota()
    {
        simd_vars _iota{};

        for (size_t offset = 0, value = 0; offset < count; ++offset)
        {
            for (size_t i = 0; i < _iota[offset].size(); ++i, ++value)
            {
                _iota[offset][i] = value;
                // std::cout << "_iota[" << offset << "][" << i << "] = i + offset * count = " << (i + offset * count) << std::endl;
            }
        }

        return _iota;
    }

    struct hasher
    {
        std::size_t operator()(simd_vars const& var) const noexcept
        {
            size_t hash;
            for (simd_int32x8 memory : var)
            {
                for (int idx : memory)
                {
                    hash = hash ^ (std::hash<int>{}(idx) << 1);
                }
            }
            return hash; // or use boost::hash_combine
        }
    };

    template <size_t count2>
    bool contains(simd_vars<count2> const & memory) const
    {
        simd_vars<count> const & memory1{*this};
        simd_vars<count2> const & memory2{memory};

        for (simd_int32x8 const & simd : memory2)
        {
            auto it = std::find(memory1.begin(), memory1.end(), simd);
            if (it == memory1.end())
                return false;
        }

        return true;
    }

    template <size_t count2>
    bool has_at_least_each_once(simd_vars<count2> const & final_memory) const
    {
        simd_vars<count> const & memory1{*this};
        simd_vars<count2> const & memory2{final_memory};

        std::unordered_map<int, int> map{};

        for (simd_int32x8 const & simd : memory1)
        {
            for (int v : simd)
                ++map[v];
        }

        for (simd_int32x8 const & simd : memory2)
        {
            for (int v : simd)
                if (!map.contains(v))
                    return false;
        }

        return true;
    }

    template <typename char_t>
    friend std::basic_ostream<char_t> & operator<<(std::basic_ostream<char_t> & cout, simd_vars<count> const & vars)
    {
        for (size_t i = 0; i < count; ++i)
        {
            simd_int32x8 const & simd = vars[i];
            cout << simd;
            if (i < count - 1)
                cout << std::endl;
        }
        return cout;
    }
};

template <size_t count>
struct memory_state
{
    void add_state(simd_vars<count> const & memory, simd_swizzle_op_trace const & trace)
    {
        auto it = _traces.find(memory);
        if (it == _traces.end())
        {
            it = _traces.insert_or_assign(it, memory, trace);
            assert(it != _traces.end());
        } else
        {
            // assert(false && "THIS SHOULD NOT HIT.");
            // merge traces
            return; // TODO: add cost model and choose smaller trace
            // std::pair<const simd_vars<1>, simd_swizzle_op_trace> & element = *it;
            // simd_swizzle_op_trace & old_traces = element.second;
            //
            // for (auto && op : trace.ops)
            // {
            //     std::visit([&](auto && op)
            //     {
            //         old_traces.ops.insert(old_traces.ops.end(), op);
            //     }, op);
            // }
        }
    }

    template <size_t count2>
    std::pair<simd_vars<count>, simd_swizzle_op_trace> find(
        simd_swizzle_state_rules const & rules,
        simd_vars<count> memory,
        simd_vars<count2> const & final_memory,
        size_t max_steps)
    {
        simd_swizzle_op_trace current_trace{};

        // add initial state
        _traces.clear();
        add_state(memory, current_trace);

        std::cout << "memory: " << memory << std::endl;
        std::cout << "_traces.size: " << _traces.size() << std::endl;

        for (size_t i = 0; i < max_steps; ++i)
        {
            if (memory.contains(final_memory))
            {
                std::cout << "search: " << final_memory << ", FOUND" << std::endl;
                break;
            }
            std::cout << "search: " << final_memory << std::endl;

            auto current_traces = _traces;
            for (auto && [old_memory, old_trace] : current_traces)
            {
                std::cout << "old_memory: " << old_memory << std::endl;
                bool found = apply_transitions(rules, old_memory, old_trace, [&](auto && new_memory, auto && new_trace)
                {
                    if (new_memory.contains(final_memory))
                    {
                        memory = new_memory;
                        current_trace = new_trace;
                        std::cout << new_memory << "; " << new_trace << ", found!!!!" << std::endl;
                        return true;
                    }
                    std::cout << new_memory << "; " << new_trace << std::endl;

                    if (new_memory.has_at_least_each_once(final_memory))
                        this->add_state(new_memory, new_trace);

                    return false;
                });

                if (found) break;
            }
        }

        return {memory, current_trace};
    }


private:
    bool apply_transitions(
        simd_swizzle_state_rules const & rules,
        simd_vars<count> const & old_memory,
        simd_swizzle_op_trace const & old_trace,
        auto && on_new_memory
    )
    {
        bool found = apply_unary_transitions(rules, old_memory, old_trace, on_new_memory);
        found = found || apply_binary_transitions(rules, old_memory, old_trace, on_new_memory);

        return found;
    }

    bool apply_unary_transitions(
        simd_swizzle_state_rules const & rules,
        simd_vars<count> const & old_memory,
        simd_swizzle_op_trace const & old_trace,
        auto && on_new_memory
    )
    {
        for (simd_swizzle_opv const & transition_rule_var : rules.transitions)
        {
            bool found = std::visit([&](auto && transition_rule)
            {
                if constexpr (transition_rule.arg_count == 1)
                {
                    for (size_t i = 0; i < count; ++i)
                    {
                        simd_int32x8 simd = old_memory[i];

                        simd = transition_rule(simd);

                        for (size_t k = 0; k < count; ++k)
                        {
                            auto current_memory = old_memory;
                            auto current_transition_rule = transition_rule;
                            auto current_trace = old_trace;

                            current_transition_rule.arg0.position = i;
                            current_transition_rule.return_argument.position = k;
                            current_trace.ops.push_back(current_transition_rule);
                            current_memory[k] = simd;
                            bool found = on_new_memory(std::move(current_memory), std::move(current_trace));
                            if (found) return true;
                        }
                    }
                }
                return false;
            }, transition_rule_var);

            if (found) return true;
        }

        return false;
    }

    bool apply_binary_transitions(
        simd_swizzle_state_rules const & rules,
        simd_vars<count> const & old_memory,
        simd_swizzle_op_trace const & old_trace,
        auto && on_new_memory
    )
    {
        for (simd_swizzle_opv const & transition_rule_var : rules.transitions)
        {
            bool found = std::visit([&](auto && transition_rule)
            {
                if constexpr (transition_rule.arg_count == 2)
                {
                    for (size_t i = 0; i < count; ++i)
                    {
                        for (size_t j = 0; j < count; ++j)
                        {
                            simd_int32x8 simd1 = old_memory[i];
                            simd_int32x8 simd2 = old_memory[j];

                            simd1 = transition_rule(simd1, simd2);

                            for (size_t k = 0; k < count; ++k)
                            {
                                auto current_memory = old_memory;
                                auto current_transition_rule = transition_rule;
                                auto current_trace = old_trace;

                                current_transition_rule.arg0.position = i;
                                current_transition_rule.arg1.position = j;
                                current_transition_rule.return_argument.position = k;
                                current_trace.ops.push_back(current_transition_rule);
                                current_memory[k] = simd1;
                                bool found = on_new_memory(std::move(current_memory), std::move(current_trace));
                                if (found) return true;
                            }
                        }
                    }
                }
                return false;
            }, transition_rule_var);

            if (found) return true;
        }

        return false;
    }

    std::unordered_map<simd_vars<count>, simd_swizzle_op_trace, typename simd_vars<count>::hasher> _traces{};
};

simd_swizzle_state_rules permute2f128_si256_transitions{
    // ------ permute2f128_si256 ------
    /*permute2f128_si256_01: A[0], A[0]*/ simd_swizzle_op{{0, 1, 2, 3,/*|*/ 0, 1, 2, 3}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_any_m256i{}, argument_imm8{0x00}},
    // /*permute2f128_si256_02: A[0], A[1]*/ simd_swizzle_op{{0, 1, 2, 3,/*|*/ 4, 5, 6, 7}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_any_m256i{}, argument_imm8{0x10}},
    /*permute2f128_si256_03: A[0], B[0]*/ simd_swizzle_op{{0, 1, 2, 3,/*|*/ 8, 9, 10, 11}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x20}},
    /*permute2f128_si256_04: A[0], B[1]*/ simd_swizzle_op{{0, 1, 2, 3,/*|*/ 12, 13, 14, 15}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x30}},
    /*permute2f128_si256_05: A[1], A[0]*/ simd_swizzle_op{{4, 5, 6, 7,/*|*/ 0, 1, 2, 3}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_any_m256i{}, argument_imm8{0x01}},
    /*permute2f128_si256_06: A[1], A[1]*/ simd_swizzle_op{{4, 5, 6, 7,/*|*/ 4, 5, 6, 7}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_any_m256i{}, argument_imm8{0x11}},
    /*permute2f128_si256_07: A[1], B[0]*/ simd_swizzle_op{{4, 5, 6, 7,/*|*/ 8, 9, 10, 11}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x21}},
    /*permute2f128_si256_08: A[1], B[1]*/ simd_swizzle_op{{4, 5, 6, 7,/*|*/ 12, 13, 14, 15}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x31}},

    // /*permute2f128_si256_09: B[0], A[0]*/ simd_swizzle_op{{8, 9, 10, 11,/*|*/ 0, 1, 2, 3}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x02}},
    // /*permute2f128_si256_10: B[0], A[1]*/ simd_swizzle_op{{8, 9, 10, 11,/*|*/ 4, 5, 6, 7}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x12}},
    // /*permute2f128_si256_11: B[0], B[0]*/ simd_swizzle_op{{8, 9, 10, 11,/*|*/ 8, 9, 10, 11}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x22}},
    // /*permute2f128_si256_12: B[0], B[1]*/ simd_swizzle_op{{8, 9, 10, 11,/*|*/ 12, 13, 14, 15}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x32}},
    // /*permute2f128_si256_13: B[1], A[0]*/ simd_swizzle_op{{12, 13, 14, 15,/*|*/ 0, 1, 2, 3}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x03}},
    // /*permute2f128_si256_14: B[1], A[1]*/ simd_swizzle_op{{12, 13, 14, 15,/*|*/ 4, 5, 6, 7}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x13}},
    // /*permute2f128_si256_15: B[1], B[0]*/ simd_swizzle_op{{12, 13, 14, 15,/*|*/ 8, 9, 10, 11}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x23}},
    // /*permute2f128_si256_16: B[1], B[1]*/ simd_swizzle_op{{12, 13, 14, 15,/*|*/ 12, 13, 14, 15}, "_mm256_permutevar8x32_epi32", argument_m256i{}, argument_m256i{}, argument_imm8{0x33}},
};

simd_swizzle_state_rules transitions{[]()
{
    simd_swizzle_state_rules _transitions{};
    auto append = [&](auto && other)
    {
        _transitions.transitions.insert(
            _transitions.transitions.end(),
            other.transitions.begin(),
            other.transitions.end()
        );
    };

    append(permute2f128_si256_transitions);

    return _transitions;
}()};

void simd_vars_test()
{
    assert((simd_vars<1>::iota() == simd_vars<1>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}}));
    assert((simd_vars<2>::iota() == simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{8+0, 8+1, 8+2, 8+3, 8+4, 8+5, 8+6, 8+7}}));
    assert((simd_vars<3>::iota() == simd_vars<3>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{8+0, 8+1, 8+2, 8+3, 8+4, 8+5, 8+6, 8+7}, simd_int32x8{16+0, 16+1, 16+2, 16+3, 16+4, 16+5, 16+6, 16+7}}));

    assert(!(simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}}.contains(simd_vars<1>{{0, 1, 2, 3, 4, 5, 7, 6}})));
    assert((simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}}.contains(simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}})));
    assert((simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{-1, -1, -1, -1, -1, -1, -1, -1}}.contains(simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}})));
    assert((simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{8, 9, 10, 11, 12, 13, 14, 15}}.contains(simd_vars<1>{{8, 9, 10, 11, 12, 13, 14, 15}})));
    assert(!(simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{8, 9, 10, 11, 12, 13, 14, 15}}.contains(simd_vars<1>{{4, 5, 6, 7, 8, 9, 10, 11}})));

    assert((simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}}.has_at_least_each_once(simd_vars<1>{{0, 1, 2, 3, 4, 5, 7, 6}})));
    assert((simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}}.has_at_least_each_once(simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}})));
    assert((simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{-1, -1, -1, -1, -1, -1, -1, -1}}.has_at_least_each_once(simd_vars<1>{{0, 1, 2, 3, 4, 5, 6, 7}})));
    assert((simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{8, 9, 10, 11, 12, 13, 14, 15}}.has_at_least_each_once(simd_vars<1>{{8, 9, 10, 11, 12, 13, 14, 15}})));
    assert((simd_vars<2>{simd_int32x8{0, 1, 2, 3, 4, 5, 6, 7}, simd_int32x8{8, 9, 10, 11, 12, 13, 14, 15}}.has_at_least_each_once(simd_vars<1>{{4, 5, 6, 7, 8, 9, 10, 11}})));
    assert(!(simd_vars<2>{simd_int32x8{0, 1, 2, 3, 5, 5, 6, 7}, simd_int32x8{8, 9, 10, 11, 12, 13, 14, 15}}.has_at_least_each_once(simd_vars<1>{{4, 5, 6, 7, 8, 9, 10, 11}})));
}

void find_identity_test()
{
    // identity
    {
        auto && [memory, final_state] = memory_state<1>{}.find(
            transitions,
            simd_vars<1>{0, 1, 2, 3, 4, 5, 6, 7}, // start state
            simd_vars<1>{0, 1, 2, 3, 4, 5, 6, 7}, // end state
            1
        );

        assert(final_state.size() == 0);
        assert((memory == simd_vars<1>{0, 1, 2, 3, 4, 5, 6, 7}));
    }
}

void find_permute2f128_si256_test()
{
    // permute2f128_si256_0x00
    {
        auto && [memory, final_state] = memory_state<1>{}.find(
            transitions,
            simd_vars<1>::iota(), // start state
            simd_vars<1>{0, 1, 2, 3, 0, 1, 2, 3}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_any_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg2.constant == 0x00);
        assert((memory == simd_vars<1>{0, 1, 2, 3, 0, 1, 2, 3}));
    }

    // permute2f128_si256_0x01
    {
        auto && [memory, final_state] = memory_state<1>{}.find(
            transitions,
            simd_vars<1>::iota(), // start state
            simd_vars<1>{4, 5, 6, 7, 0, 1, 2, 3}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_any_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg2.constant == 0x01);
        assert((memory == simd_vars<1>{4, 5, 6, 7, 0, 1, 2, 3}));
    }

    // permute2f128_si256_0x11
    {
        auto && [memory, final_state] = memory_state<1>{}.find(
            transitions,
            simd_vars<1>::iota(), // start state
            simd_vars<1>{4, 5, 6, 7, 4, 5, 6, 7}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_any_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg2.constant == 0x11);
        assert((memory == simd_vars<1>{4, 5, 6, 7, 4, 5, 6, 7}));
    }

    // permute2f128_si256_0x20
    {
        std::cout << "~~~~~~~~~~" << std::endl;
        auto && [memory, final_state] = memory_state<2>{}.find(
            transitions,
            simd_vars<2>::iota(), // start state
            simd_vars<1>{0, 1, 2, 3, 8, 9, 10, 11}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg1.position == 1);
        assert(state0.arg2.constant == 0x20);
        assert((memory == simd_vars<2>{simd_int32x8{0,1,2,3,8,9,10,11}, simd_int32x8{8,9,10,11,12,13,14,15}}));
    }

    // permute2f128_si256_0x30
    {
        std::cout << "~~~~~~~~~~" << std::endl;
        auto && [memory, final_state] = memory_state<2>{}.find(
            transitions,
            simd_vars<2>::iota(), // start state
            simd_vars<1>{0, 1, 2, 3, 12, 13, 14, 15}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg1.position == 1);
        assert(state0.arg2.constant == 0x30);
        assert((memory == simd_vars<2>{simd_int32x8{0, 1, 2, 3, 12, 13, 14, 15}, simd_int32x8{8,9,10,11,12,13,14,15}}));
    }

    // permute2f128_si256_0x31
    {
        std::cout << "~~~~~~~~~~" << std::endl;
        auto && [memory, final_state] = memory_state<2>{}.find(
            transitions,
            simd_vars<2>::iota(), // start state
            simd_vars<1>{4, 5, 6, 7, 12, 13, 14, 15}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 0);
        assert(state0.arg1.position == 1);
        assert(state0.arg2.constant == 0x31);
        assert((memory == simd_vars<2>{simd_int32x8{4, 5, 6, 7, 12, 13, 14, 15}, simd_int32x8{8,9,10,11,12,13,14,15}}));
    }

    // permute2f128_si256_0x20 reverse
    {
        std::cout << "~~~~~~~~~~" << std::endl;
        auto && [memory, final_state] = memory_state<2>{}.find(
            transitions,
            simd_vars<2>::iota(), // start state
            simd_vars<1>{8, 9, 10, 11, 0, 1, 2, 3}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 1);
        assert(state0.arg1.position == 0);
        assert(state0.arg2.constant == 0x20);
        assert((memory == simd_vars<2>{simd_int32x8{8, 9, 10, 11, 0, 1, 2, 3}, simd_int32x8{8,9,10,11,12,13,14,15}}));
    }

    // permute2f128_si256_0x21 reverse
    {
        std::cout << "~~~~~~~~~~" << std::endl;
        auto && [memory, final_state] = memory_state<2>{}.find(
            transitions,
            simd_vars<2>::iota(), // start state
            simd_vars<1>{12, 13, 14, 15, 0, 1, 2, 3}, // end state
            1
        );

        assert(final_state.size() == 1);
        auto state0 = std::get<simd_swizzle_op<argument_m256i, argument_m256i, argument_imm8>>(final_state.ops[0]);
        assert(state0.name == "_mm256_permutevar8x32_epi32");
        assert(state0.return_argument.position == 0);
        assert(state0.arg0.position == 1);
        assert(state0.arg1.position == 0);
        assert(state0.arg2.constant == 0x21);
        assert((memory == simd_vars<2>{simd_int32x8{12, 13, 14, 15, 0, 1, 2, 3}, simd_int32x8{8,9,10,11,12,13,14,15}}));
    }
}

int main()
{
    simd_vars_test();
    find_identity_test();
    find_permute2f128_si256_test();
}
