
#pragma once

#include <cmath>
#include <limits>
#include <iostream>
#include <vector>

namespace stellar::test
{

template <auto compute_forward_fn, auto compute_backward_fn>
void compute_minimiser_blocks_test()
{
    std::vector<int> values{
        3 ,8 ,5 ,4 ,10,17,
        9 ,15,20,13,19,13,
        0 ,1 ,4 ,0 ,8 ,18,
        8 ,16,14,9 ,14,8 ,
        19,3 ,3 ,9 ,5 ,20,
        0 ,1 ,2 ,3 ,4 ,5 ,
        6 ,7 ,8 ,9 ,10,11,
        11,10,9 ,8 ,7 ,6 ,
        // 5, 4, 3, 2, 1, 0
    };

    std::vector<int> window_value_forward_expected{
        /*1:*/ 3 /*3 */,4 /*8 */,4 /*5 */,4 /*4 */,10/*10*/,
        /*2:*/ 9 /*17*/,9 /*9 */,13/*15*/,13/*20*/,13/*13*/,
        /*3:*/ 0 /*19*/,0 /*13*/,0 /*0 */,1 /*1 */,4 /*4 */,
        /*4:*/ 0 /*0 */,8 /*8 */,8 /*18*/,8 /*8 */,16/*16*/,
        /*5:*/ 8 /*14*/,8 /*9 */,8 /*14*/,8 /*8 */,19/*19*/,
        /*6:*/ 3 /*3 */,3 /*3 */,5 /*9 */,5 /*5 */,20/*20*/,
        /*7:*/ 0 /*0 */,1 /*1 */,2 /*2 */,3 /*3 */,4 /*4 */,
        /*8:*/ 5 /*5 */,6 /*6 */,7 /*7 */,8 /*8 */,9 /*9 */,
        /*9:*/ 9 /*10*/,9 /*11*/,9 /*11*/,9 /*10*/,9 /*9 */,
        /*0:*/ 6 /*8 */,6 /*7 */,6 /*6 */,
    };
    std::vector<int> window_offset_forward_expected{
        /*1:*/ 0, 3, 3, 3, 4,
        /*2:*/ 1, 1, 4, 4, 4,
        /*3:*/ 2, 2, 2, 3, 4,
        /*4:*/ 0, 3, 3, 3, 4,
        /*5:*/ 3, 3, 3, 3, 4,
        /*6:*/ 1, 1, 3, 3, 4,
        /*7:*/ 0, 1, 2, 3, 4,
        /*8:*/ 0, 1, 2, 3, 4,
        /*9:*/ 4, 4, 4, 4, 4,
        /*0:*/ 2, 2, 2,
    };

    std::vector<int> window_value_backward_expected{
        /*1:*/ 3 /*3 */,3 /*8 */,3 /*5 */,3 /*4 */,3 /*10*/,
        /*2:*/ 17/*17*/,9 /*9 */,9 /*15*/,9 /*20*/,9 /*13*/,
        /*3:*/ 19/*19*/,13/*13*/,0 /*0 */,0 /*1 */,0 /*4 */,
        /*4:*/ 0 /*0 */,0 /*8 */,0 /*18*/,0 /*8 */,0 /*16*/,
        /*5:*/ 14/*14*/,9 /*9 */,9 /*14*/,8 /*8 */,8 /*19*/,
        /*6:*/ 3 /*3 */,3 /*3 */,3 /*9 */,3 /*5 */,3 /*20*/,
        /*7:*/ 0 /*0 */,0 /*1 */,0 /*2 */,0 /*3 */,0 /*4 */,
        /*8:*/ 5 /*5 */,5 /*6 */,5 /*7 */,5 /*8 */,5 /*9 */,
        /*9:*/ 10/*10*/,10/*11*/,10/*11*/,10/*10*/,9 /*9 */,
        /*0:*/ 8 /*8 */,7 /*7 */,6 /*6 */,
    };
    std::vector<int> window_offset_backward_expected{
        /*1:*/ 0, 0, 0, 0, 0,
        /*2:*/ 0, 1, 1, 1, 1,
        /*3:*/ 0, 1, 2, 2, 2,
        /*4:*/ 0, 0, 0, 0, 0,
        /*5:*/ 0, 1, 1, 3, 3,
        /*6:*/ 0, 1, 1, 1, 1,
        /*7:*/ 0, 0, 0, 0, 0,
        /*8:*/ 0, 0, 0, 0, 0,
        /*9:*/ 0, 0, 0, 3, 4,
        /*0:*/ 0, 1, 2,
        // 5, 5, 5, 5, 5, 5,
    };

    std::size_t window_size = 6;

    std::vector<int> window_value_forward;
    std::vector<int> window_offset_forward;
    window_value_forward.resize(values.size());
    window_offset_forward.resize(values.size());

    std::vector<int> window_value_backward;
    std::vector<int> window_offset_backward;
    window_value_backward.resize(values.size());
    window_offset_backward.resize(values.size());

    compute_forward_fn(window_size, values.data(), window_value_forward.data(), window_value_forward.size(), window_offset_forward.data());
    compute_backward_fn(window_size, values.data(), window_value_backward.data(), window_value_backward.size(), window_offset_backward.data());

    std::cout << "window_value_forward: " << std::endl;
    print_chunked(window_value_forward, window_size - 1);
    assert(std::equal(window_value_forward_expected.begin(), window_value_forward_expected.end(), window_value_forward.begin()));

    std::cout << "window_offset_forward: " << std::endl;
    print_chunked(window_offset_forward, window_size - 1);
    assert(std::equal(window_offset_forward_expected.begin(), window_offset_forward_expected.end(), window_offset_forward.begin()));

    std::cout << "window_value_backward: " << std::endl;
    print_chunked(window_value_backward, window_size - 1);
    assert(std::equal(window_value_backward_expected.begin(), window_value_backward_expected.end(), window_value_backward.begin()));

    std::cout << "window_offset_backward: " << std::endl;
    print_chunked(window_offset_backward, window_size - 1);
    assert(std::equal(window_offset_backward_expected.begin(), window_offset_backward_expected.end(), window_offset_backward.begin()));
}

} // namespace stellar::test
