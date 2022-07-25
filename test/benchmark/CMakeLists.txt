
include (seqan3_require_benchmark)
seqan3_require_benchmark()

macro (add_benchmark test_filename)
    seqan3_test_component (target "${test_filename}" TARGET_NAME)

    add_executable("${target}" "${test_filename}")

    target_link_libraries ("${target}" PUBLIC benchmark_main)
    target_link_libraries ("${target}" PUBLIC "${PROJECT_NAME}_lib")
    # add_dependencies (benchmark ${target})
endmacro ()

add_subdirectory (prefilter/local_query_prefilter/minimiser_window)
