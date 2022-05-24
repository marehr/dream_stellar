# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/stellar/CMakeLists.txt
#
# CMakeLists.txt file for STELLAR.
# ===========================================================================

cmake_minimum_required (VERSION 3.10)
message (STATUS "Configuring apps/stellar")

set (SEQAN_APP_VERSION "${PROJECT_VERSION}")

# Compile times with `make -j 6` measured on gcc 11.1:
#
# `cmake -DSTELLAR_PARALLEL_BUILD=OFF <...>`: (default behaviour!)
#     is recommend for a single time build. It will produce the fastest and smallest
#     binary of all available options and is therefore the **default**.
#
# `cmake -DSTELLAR_PARALLEL_BUILD=ON <...>`:
#     will increase runtime by around 4-5% in RELEASE mode,
#     but will be 30% (DEBUG) - 45% (RELEASE) faster on parallel rebuilds.
#     For a single RELEASE build it can speed up the compilation by 22%.
#
# `cmake -DSTELLAR_PARALLEL_BUILD=ON -DSTELLAR_PRECOMPILED_HEADER_BUILD=ON <...>`:
#     will increase runtime by around 4-5% in RELEASE mode,
#     but will be 25% (DEBUG) - 48% (RELEASE) faster on parallel rebuilds.
#     For a single RELEASE build it can speed up the compilation by 21%.
#
# #### Cold build ####
#
# i.e. everything must be build:
#
# | CMAKE_BUILD_TYPE | -DSTELLAR_PARALLEL_BUILD=        | wall time | usr time | system | cpu% |
# | ---------------- | -------------------------------- | --------- | -------- | ------ | ---- |
# | DEBUG            | OFF                              | 28.22s    | 27.10s   | 1.02s  |  99% |
# | DEBUG            | ON                               | 30.97s    | 115.83s  | 4.08s  | 387% |
# | DEBUG            | ON & PRECOMPILED_HEADER_BUILD=ON | 30.88s    | 78.85s   | 3.80s  | 267% |
# | RELEASE          | OFF                              | 51.68s    | 50.72s   | 0.82s  |  99% |
# | RELEASE          | ON                               | 40.69s    | 152.75s  | 3.24s  | 383% |
# | RELEASE          | ON & PRECOMPILED_HEADER_BUILD=ON | 40.85s    | 121.12s  | 3.08s  | 304% |
#
# #### Partial rebuild ####
#
# i.e. only stellar.main.tpp is updated, everything is left untouched.
#
# | CMAKE_BUILD_TYPE | -DSTELLAR_PARALLEL_BUILD=        | wall time | usr time | system | cpu% |
# | ---------------- | -------------------------------- | --------- | -------- | ------ | ---- |
# | DEBUG            | OFF                              | 28.22s    | 27.10s   | 1.02s  |  99% | (same as before)
# | DEBUG            | ON                               | 20.51s    | 100.46s  | 3.43s  | 506% |
# | DEBUG            | ON & PRECOMPILED_HEADER_BUILD=ON | 21.39s    | 70.68s   | 3.43s  | 346% |
# | RELEASE          | OFF                              | 51.68s    | 50.72s   | 0.82s  |  99% | (same as before)
# | RELEASE          | ON                               | 27.41s    | 146.22s  | 3.05s  | 544% |
# | RELEASE          | ON & PRECOMPILED_HEADER_BUILD=ON | 26.94s    | 113.51s  | 2.68s  | 431% |
option (STELLAR_PARALLEL_BUILD "Compile stellar on multiple translation units. (Enables faster build when building in parallel `make -j 6`)" OFF)

include (CMakeDependentOption)
cmake_dependent_option (STELLAR_PRECOMPILED_HEADER_BUILD "Use pre-compiled header to speed up (single) build time." On "STELLAR_PARALLEL_BUILD" OFF)

# ----------------------------------------------------------------------------
# Setup SeqAn 2
# ----------------------------------------------------------------------------

find_path (STELLAR_ROOT_DIR NAMES src/stellar.cpp HINTS "${CMAKE_CURRENT_LIST_DIR}/..")

# needed for find_package(SeqAn) (finding header files)
find_path (SEQAN_INCLUDE_PATH NAMES seqan/version.h HINTS "${STELLAR_ROOT_DIR}/lib/seqan/include")

find_path (SEQAN_BUILD_SYSTEM_CMAKE_DIR NAMES SeqAnBuildSystem.cmake HINTS "${STELLAR_ROOT_DIR}/lib/seqan/util/cmake/")
if (EXISTS "${SEQAN_BUILD_SYSTEM_CMAKE_DIR}")
    # needed for include(SeqAnBuildSystem)
    set (CMAKE_MODULE_PATH "${SEQAN_BUILD_SYSTEM_CMAKE_DIR}" "${CMAKE_MODULE_PATH}")
    # needed for find_package (SeqAn) (finding seqan-config.cmake)
    set (CMAKE_PREFIX_PATH "${SEQAN_BUILD_SYSTEM_CMAKE_DIR}" "${CMAKE_PREFIX_PATH}")
endif ()

add_custom_target (seqan_library) # place-holder needed for SeqAnBuildSystem
include (SeqAnBuildSystem) # for seqan_setup_install_vars

# ----------------------------------------------------------------------------
# Setup SeqAn 3
# ----------------------------------------------------------------------------

# A interface to reuse common properties.
# You can add more external include paths of other projects that are needed for your project.
add_library ("${PROJECT_NAME}_interface" INTERFACE)
# target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE seqan3::seqan3)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE ../include)
target_compile_options ("${PROJECT_NAME}_interface" INTERFACE "-pedantic" "-Wall" "-Wextra")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (SeqAn CONFIG REQUIRED)
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
target_include_directories("${PROJECT_NAME}_interface" INTERFACE ${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
target_compile_definitions ("${PROJECT_NAME}_interface" INTERFACE ${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
target_compile_definitions ("${PROJECT_NAME}_interface" INTERFACE -DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
target_compile_definitions ("${PROJECT_NAME}_interface" INTERFACE -DSEQAN_REVISION="${SEQAN_REVISION}")
target_compile_definitions ("${PROJECT_NAME}_interface" INTERFACE -DSEQAN_DATE="${SEQAN_DATE}")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Require at least C++-17
target_compile_features ("${PROJECT_NAME}_interface" INTERFACE cxx_std_17)

# Require OpenMP
find_package (OpenMP) # requires cmake 3.10
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE OpenMP::OpenMP_CXX)

if (STELLAR_PARALLEL_BUILD)
    target_compile_definitions ("${PROJECT_NAME}_interface" INTERFACE -DSTELLAR_PARALLEL_BUILD)

    add_library (stellar_main OBJECT
                 stellar/stellar.main.dna.cpp
                 stellar/stellar.main.dna5.cpp
                 stellar/stellar.main.rna.cpp
                 stellar/stellar.main.rna5.cpp
                 stellar/stellar.main.aminoacid.cpp
                 stellar/stellar.main.char.cpp)
    target_link_libraries (stellar_main PUBLIC "${PROJECT_NAME}_interface")

    add_library (stellar_arg_parser OBJECT stellar/stellar.arg_parser.cpp)
    target_link_libraries (stellar_arg_parser PUBLIC "${PROJECT_NAME}_interface")

    add_library (stellar_diagnostics OBJECT stellar/stellar.diagnostics.cpp)
    target_link_libraries (stellar_diagnostics PUBLIC "${PROJECT_NAME}_interface")
endif ()

if (STELLAR_PARALLEL_BUILD AND STELLAR_PRECOMPILED_HEADER_BUILD)
    target_precompile_headers (stellar_main PRIVATE ../include/stellar/app/stellar.main.tpp)
endif ()

# Update the list of file names below if you add source files to your application.
add_executable (stellar stellar.cpp)

if (STELLAR_PARALLEL_BUILD)
    target_link_libraries (stellar PRIVATE stellar_main stellar_arg_parser stellar_diagnostics)
endif ()

# Add dependencies found by find_package (SeqAn).
target_link_libraries (stellar PUBLIC "${PROJECT_NAME}_interface" ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (stellar)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install stellar in ${PREFIX}/bin directory
install (TARGETS stellar
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/stellar for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/NC_001474.fasta
               example/NC_001477.fasta
               example/reads.fasta
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:stellar")
  set (CPACK_PACKAGE_NAME "stellar")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "STELLAR - SwifT Exact LocaL AligneR")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (stellar "Stellar")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:stellar")
