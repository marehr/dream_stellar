#include <cassert>
#include <vector>

#include "simd_fn.hpp"

void gather_test()
{
    size_t window_size = 6;
    std::vector<int> values{
        3, 8, 5, 4, 10,17,
        9, 15,20,13,19,13,
        0, 1, 4, 0, 8, 18,
        8, 16,14,9, 14,8,
        19,3, 3, 9, 5, 20,
        0, 1, 2, 3, 4, 5,
        6, 7, 8, 9, 10,11,
        11,10,9, 8, 7, 6};

    {
        int32x8_t column1 = simd_gather(window_size, 0, values.data());
        int column1_expected[]{3, 9, 0, 8, 19, 0, 6, 11};
        assert(std::equal(column1_expected, column1_expected + 8, (int *)(&column1[0])));
    }

    {
        int32x8_t column2 = simd_gather(window_size, 1, values.data());
        int column2_expected[]{8, 15, 1, 16, 3, 1, 7, 10};
        assert(std::equal(column2_expected, column2_expected + 8, (int *)(&column2[0])));
    }

    {
        int32x8_t column3 = simd_gather(window_size, 2, values.data());
        int column3_expected[]{5, 20, 4, 14, 3, 2, 8, 9};
        assert(std::equal(column3_expected, column3_expected + 8, (int *)(&column3[0])));
    }

    {
        int32x8_t column4 = simd_gather(window_size, 3, values.data());
        int column4_expected[]{4, 13, 0, 9, 9, 3, 9, 8};
        assert(std::equal(column4_expected, column4_expected + 8, (int *)(&column4[0])));
    }

    {
        int32x8_t column5 = simd_gather(window_size, 4, values.data());
        int column5_expected[]{10, 19, 8, 14, 5, 4, 10, 7};
        assert(std::equal(column5_expected, column5_expected + 8, (int *)(&column5[0])));
    }

    {
        int32x8_t column6 = simd_gather(window_size, 5, values.data());
        int column6_expected[]{17, 13, 18, 8, 20, 5, 11, 6};
        assert(std::equal(column6_expected, column6_expected + 8, (int *)(&column6[0])));
    }
}

int main()
{
    gather_test();

    return 0;
}
