
add_library (stellar_local_query_prefilter INTERFACE)

set (DREAM_STELLAR_BUILD_THRESHOLD 1)
if (DREAM_STELLAR_BUILD_THRESHOLD)
    find_path (RAPTOR_THRESHOLD_INCLUDE_DIR NAMES raptor/threshold/threshold.hpp HINTS "../lib/valik/lib/raptor/include")
    find_path (RAPTOR_THRESHOLD_SRC NAMES threshold/precompute_threshold.cpp HINTS "../lib/valik/lib/raptor/src")
    find_path (VALIK_INCLUDE_DIR NAMES valik/valik.hpp HINTS "../lib/valik/include")

    target_include_directories ("${PROJECT_NAME}_interface" INTERFACE )

    add_library (raptor_interface INTERFACE)
    target_include_directories (raptor_interface INTERFACE "${RAPTOR_THRESHOLD_INCLUDE_DIR}")
    target_link_libraries (raptor_interface INTERFACE seqan3::seqan3)

    add_subdirectory ("${RAPTOR_THRESHOLD_SRC}/threshold" "raptor_threshold_build")

    #  "stellar/prefilter/local_query_prefilter.cpp"
    target_include_directories (stellar_local_query_prefilter INTERFACE "${VALIK_INCLUDE_DIR}")
    target_link_libraries (stellar_local_query_prefilter INTERFACE raptor_threshold)
endif ()
