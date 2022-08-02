#include <cassert>
#include <iostream>
#include <limits>
#include <span>
#include <vector>

#include "simd_minimiser_window.hpp"

#include "diagnostics.hpp"
#include "compute_minimiser_blocks.hpp"
#include "compute_minimiser_blocks_simd.hpp"
#include "compute_minimiser_blocks_test.hpp"

using namespace stellar::test;

int main()
{
    gather_test();
    compute_minimiser_blocks_test<compute_forward_full, compute_backward_full>();
    minimiser_test();
}
