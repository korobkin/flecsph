# Macro to add googletest
macro(package_add_test TESTNAME)
    add_executable(${TESTNAME} ${ARGN} ${FleCSI_RUNTIME}/runtime_driver.cc)
    target_include_directories(${TESTNAME}
        PRIVATE
            ${CMAKE_SOURCE_DIR}/mpisph
            ${CMAKE_SOURCE_DIR}/app/drivers/include
    )

    target_link_libraries(${TESTNAME}
        PRIVATE
            flecsph::library_flags
            flecsph::compile_flags
            FleCSI::flecsi
    )
    add_test(
        NAME ${TESTNAME}
        COMMAND ${TESTNAME}
        WORKING_DIRECTORY ${PROJECT_DIR}
    )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

# Macro to add googletest for MPI
macro(package_add_test_MPI TESTNAME)
    add_executable(${TESTNAME} ${ARGN} ${FleCSI_RUNTIME}/runtime_driver.cc)
    target_include_directories(${TESTNAME}
        PRIVATE
            ${CMAKE_SOURCE_DIR}/mpisph
            ${CMAKE_SOURCE_DIR}/app/drivers/include
    )
    target_link_libraries(${TESTNAME}
        PRIVATE
            flecsph::library_flags
            flecsph::compile_flags
            FleCSI::flecsi
    )
    add_test(
        NAME ${TESTNAME}
        COMMAND ${MPIEXEC} -n 4 ${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME}
        WORKING_DIRECTORY ${PROJECT_DIR}
    )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()