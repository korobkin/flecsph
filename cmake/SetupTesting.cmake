# Macro to add googletest
macro(package_add_test TESTNAME PARNAME DRIVER)
    add_test(
        NAME ${TESTNAME}
        COMMAND ${MPIEXEC} -n 1 ${PROJECT_BINARY_DIR}/app/drivers/${DRIVER} ${PROJECT_SOURCE_DIR}/data/${PARNAME}
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/tests"
    )
endmacro()

# Macro to add googletest for MPI
macro(package_add_test_MPI TESTNAME PARNAME DRIVER)
    add_test(
        NAME ${TESTNAME}
        COMMAND ${MPIEXEC} -n 4 ${PROJECT_BINARY_DIR}/app/drivers/${DRIVER} ${PROJECT_SOURCE_DIR}/data/${PARNAME}
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/tests"
    )
endmacro()

# Macro to add googletest
macro(package_add_test_executable TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    target_link_libraries(${TESTNAME}
        PRIVATE
            flecsph::flags
    )
    add_test(
        NAME ${TESTNAME}
        #COMMAND ${TESTNAME}
        COMMAND ${MPIEXEC} -n 1 "${PROJECT_BINARY_DIR}/tests/${TESTNAME}"
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/tests"
    )
    #set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    set_target_properties(${TESTNAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/tests")
endmacro()

# Macro to add googletest for MPI
macro(package_add_test_executable_MPI TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    target_link_libraries(${TESTNAME}
        PRIVATE
            flecsph::flags
    )
    add_test(
        NAME ${TESTNAME}
        COMMAND ${MPIEXEC} -n 4 "${PROJECT_BINARY_DIR}/tests/${TESTNAME}"
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/tests"
    )
    #set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    set_target_properties(${TESTNAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/tests")
endmacro()
