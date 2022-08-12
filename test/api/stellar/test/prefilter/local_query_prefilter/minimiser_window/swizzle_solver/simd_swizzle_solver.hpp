
#pragma once

#include <chrono>
#include <unordered_map>

#include "simd_swizzle_op_trace.hpp"
#include "simd_swizzle_state_rules.hpp"
#include "simd_vars.hpp"

template <size_t count>
struct simd_swizzle_solver
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
    std::pair<simd_vars<count>, simd_swizzle_op_trace> solve(
        simd_swizzle_state_rules const & rules,
        simd_vars<count> memory,
        simd_vars<count2> const & final_memory,
        size_t max_steps)
    {
        simd_swizzle_op_trace current_trace{};

        // add initial state
        _traces.clear();
        _traces.reserve(5 * 1024 * 1024 * 1024 / sizeof(simd_swizzle_opv));
        add_state(memory, current_trace);

        // std::cout << "memory: " << memory << std::endl;
        // std::cout << "_traces.size: " << _traces.size() << std::endl;

        for (size_t i = 0; i < max_steps; ++i)
        {
            if (memory.contains(final_memory))
            {
                std::cout << "search: " << final_memory << ", FOUND" << std::endl;
                break;
            }
            // std::cout << "search: " << final_memory << std::endl;

            std::cout << i << "): |traces| = " << _traces.size() << std::endl;

            auto start = std::chrono::steady_clock::now();
            size_t k = 0;
            auto current_traces = _traces;
            for (auto && [old_memory, old_trace] : current_traces)
            {
                size_t progress = ((k * 100) / (current_traces.size()));
                size_t progress_tick = ((k * 100) % (current_traces.size()));
                if (progress_tick < 100){
                    auto end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> elapsed_seconds = end-start;
                    start = end;
                    std::cout << i << ") " << progress << "/ 100: " << elapsed_seconds.count() << "s" << std::endl;
                }
                ++k;

                // std::cout << "old_memory: " << old_memory << std::endl;
                bool found = apply_transitions(rules, old_memory, old_trace, [&](auto && new_memory, auto && new_trace)
                {
                    if (new_memory.contains(final_memory))
                    {
                        memory = new_memory;
                        current_trace = new_trace;
                        std::cout << new_memory << "; " << new_trace << ", found!!!!" << std::endl;
                        return true;
                    }
                    // std::cout << new_memory << "; " << new_trace << std::endl;

                    if (new_memory.has_at_least_each_once(final_memory))
                        this->add_state(new_memory, new_trace);

                    return false;
                });

                if (found) break;
            }
            std::cout << std::endl;
        }

        return {memory, current_trace};
    }


private:
    static bool apply_transitions(
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

    static bool apply_unary_transitions(
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

    static bool apply_binary_transitions(
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
