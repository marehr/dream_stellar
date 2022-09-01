#include "transpose_matrix.hpp"

using namespace stellar::test;

int main()
{
    std::cout << "transpose_matrix_int32x4x4_test<transpose_matrix_32x4x4_sse4>" << std::endl;
    transpose_matrix_int32x4x4_test<transpose_matrix_32x4x4_sse4>();
    std::cout << "transpose_matrix_int32x4x4_test<transpose_matrix_32x4x4_avx2>" << std::endl;
    transpose_matrix_int32x4x4_test<transpose_matrix_32x4x4_avx2>();
    std::cout << "transpose_matrix_int32x4x4_test<transpose_matrix_32x4x4_avx2_gather>" << std::endl;
    transpose_matrix_int32x4x4_test<transpose_matrix_32x4x4_avx2_gather>();
    std::cout << "transpose_matrix_int32x4x4_test<transpose_matrix_32x4x4_omp>" << std::endl;
    transpose_matrix_int32x4x4_test<transpose_matrix_32x4x4_omp>();
    std::cout << std::endl;

    std::cout << "transpose_matrix_int32x8x8_test<transpose_matrix_32x8x8_avx2>" << std::endl;
    transpose_matrix_int32x8x8_test<transpose_matrix_32x8x8_avx2>();
    std::cout << "transpose_matrix_int32x8x8_test<transpose_matrix_32x8x8_avx2_gather>" << std::endl;
    transpose_matrix_int32x8x8_test<transpose_matrix_32x8x8_avx2_gather>();
    std::cout << "transpose_matrix_int32x8x8_test<transpose_matrix_32x8x8_omp>" << std::endl;
    transpose_matrix_int32x8x8_test<transpose_matrix_32x8x8_omp>();
    std::cout << std::endl;

    std::cout << "transpose_matrix_int32x6x4_test<transpose_matrix_32x6x4_omp>" << std::endl;
    transpose_matrix_int32x6x4_test<transpose_matrix_32x6x4_omp>();
    std::cout << "transpose_matrix_int32x6x4_test<transpose_matrix_32x6x4_sse4>" << std::endl;
    transpose_matrix_int32x6x4_test<transpose_matrix_32x6x4_sse4>();
    std::cout << "transpose_matrix_int32x6x4_test<transpose_matrix_32x6x4_avx2>" << std::endl;
    transpose_matrix_int32x6x4_test<transpose_matrix_32x6x4_avx2>();
    std::cout << "transpose_matrix_int32x6x4_test<transpose_matrix_32x6x4_avx2_gather>" << std::endl;
    transpose_matrix_int32x6x4_test<transpose_matrix_32x6x4_avx2_gather>();
    std::cout << std::endl;
}
