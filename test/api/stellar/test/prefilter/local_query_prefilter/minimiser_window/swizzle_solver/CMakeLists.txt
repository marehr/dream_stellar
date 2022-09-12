
add_executable (swizzle_solver swizzle_solver.cpp)
target_include_directories (swizzle_solver PUBLIC ../)
target_compile_features (swizzle_solver PUBLIC cxx_std_20)
