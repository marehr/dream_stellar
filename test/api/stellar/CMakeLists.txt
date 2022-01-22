cmake_minimum_required (VERSION 3.8)

add_subdirectory (prefilter)
add_subdirectory (test)

add_api_test (seqan2to3_alphabet_test.cpp)
add_api_test (stellar_database_segment_test.cpp)
add_api_test (stellar_index_test.cpp)
add_api_test (query_matches_test.cpp)
