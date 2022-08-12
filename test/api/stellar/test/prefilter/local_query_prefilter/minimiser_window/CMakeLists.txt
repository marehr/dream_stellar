
add_api_test (minimiser_window_test.cpp)
add_api_test (simd_minimiser_window_test_detail.cpp)

target_compile_options(simd_minimiser_window_test_detail PUBLIC "-march=skylake")

add_api_test (swizzle_solver.cpp)
