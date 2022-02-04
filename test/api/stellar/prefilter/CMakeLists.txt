cmake_minimum_required (VERSION 3.8)

add_subdirectory (detail)

add_api_test (database_agent_distributer_test.cpp)
add_api_test (no_prefilter_test.cpp)
add_api_test (nsegment_database_agent_splitter_test.cpp)
add_api_test (prefilter_agent_test.cpp)
add_api_test (prefilter_agents_test.cpp)
add_api_test (whole_database_agent_splitter_test.cpp)
