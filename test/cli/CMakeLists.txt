cmake_minimum_required (VERSION 3.8)

add_custom_target (cli_test ALL)
add_dependencies (cli_test stellar)

# ----------------------------------------------------------------------------
# Setup SeqAn 2 Tests
# ----------------------------------------------------------------------------
find_path (SEQAN_APP_TEST_PYTHON_LIB_DIR NAMES seqan/app_tests.py HINTS "${STELLAR_ROOT_DIR}/lib/seqan/util/py_lib/")
macro (stellar_cli_test APP_NAME alphabet database_strand output_extension)
    find_package (PythonInterp 3)
    if (PYTHONINTERP_FOUND)
      add_test (NAME "cli/${APP_NAME}/${alphabet}/${database_strand}/${output_extension}"
                COMMAND ${CMAKE_COMMAND} -E env
                            PYTHONPATH=${SEQAN_APP_TEST_PYTHON_LIB_DIR} # make seqan/app_tests.py available
                        ${PYTHON_EXECUTABLE}
                        ${STELLAR_ROOT_DIR}/test/cli/run_tests.py
                        ${STELLAR_ROOT_DIR} ${CMAKE_BINARY_DIR}
                        --alphabets "${alphabet}"
                        --database-strands "${database_strand}"
                        --output-extensions "${output_extension}"
                )
    endif (PYTHONINTERP_FOUND)
endmacro (stellar_cli_test APP_NAME)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

foreach(alphabet dna dna5 protein char)
    foreach(database_strand forward reverse both)
        foreach(output_extensions gff txt)
            stellar_cli_test (stellar "${alphabet}" "${database_strand}" "${output_extensions}")
        endforeach()
    endforeach()
endforeach()


# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable stellar in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} stellar CACHE INTERNAL "")

#Add test to check correct ctd generation.
seqan_add_ctd_test (stellar)
