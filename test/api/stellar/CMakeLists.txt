cmake_minimum_required (VERSION 3.8)

add_subdirectory (extension)
add_subdirectory (prefilter)
add_subdirectory (test)
add_subdirectory (utils)
add_subdirectory (verification)

add_api_test (stellar_database_segment_test.cpp)
add_api_test (stellar_index_test.cpp)
