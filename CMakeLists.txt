#------------------------------------------------------------------------------#
# Copyright (c) 2020 Triad National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

project(flecsph)
cmake_minimum_required(VERSION 3.12)

#------------------------------------------------------------------------------#
# The project name
#------------------------------------------------------------------------------#
project(flecsph)

# Macro to add googletest 
macro(package_add_test TESTNAME)
    add_executable(${TESTNAME} ${ARGN} ${FleCSI_RUNTIME}/runtime_driver.cc)
    target_link_libraries(${TESTNAME} ${FleCSPH_LIBRARIES})
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
    target_link_libraries(${TESTNAME} ${FleCSPH_LIBRARIES})
    add_test(
        NAME ${TESTNAME}
        COMMAND ${MPIEXEC} -n 4 ${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME}
        WORKING_DIRECTORY ${PROJECT_DIR}
    )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

#------------------------------------------------------------------------------#
# We require a C++17-capable compiler
#------------------------------------------------------------------------------#
if(CYGWIN AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif()

#------------------------------------------------------------------------------#
# Find FleCSI
#------------------------------------------------------------------------------#
find_package(FleCSI REQUIRED)
include_directories(${FleCSI_INCLUDE_DIR})
message(STATUS ${FleCSI_INCLUDE_DIR})
set(FleCSPH_LIBRARIES ${FleCSI_LIBRARIES})
message(STATUS ${FleCSPH_LIBRARIES})

#------------------------------------------------------------------------------#
# Find GoogleTest
#------------------------------------------------------------------------------#
enable_testing()
find_package(GTest REQUIRED)
list (APPEND FleCSPH_LIBRARIES GTest::GTest GTest::Main)

#------------------------------------------------------------------------------#
# Find OpenMP 
#------------------------------------------------------------------------------#
find_package(OpenMP REQUIRED)
list (APPEND FleCSPH_LIBRARIES OpenMP::OpenMP_CXX)

#------------------------------------------------------------------------------#
# MPI
#------------------------------------------------------------------------------#
find_package(MPI REQUIRED)
list (APPEND FleCSPH_LIBRARIES MPI::MPI_CXX )

#------------------------------------------------------------------------------#
# Add GSL
#------------------------------------------------------------------------------#
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
list(APPEND FleCSPH_LIBRARIES ${GSL_LIBRARIES})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GSL_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GSL_CXX_FLAGS}")

#------------------------------------------------------------------------------#
# Add Boost
#------------------------------------------------------------------------------#
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

#------------------------------------------------------------------------------#
# HDF5
#------------------------------------------------------------------------------#
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED)
list(APPEND FleCSPH_LIBRARIES ${HDF5_LIBRARIES})
include_directories(${HDF5_INCLUDE_DIR})

#------------------------------------------------------------------------------#
# Add mpisph tests
#------------------------------------------------------------------------------#
add_subdirectory("mpisph/")
add_subdirectory("include/physics/test")
add_subdirectory("include/tree_topology/test")

#------------------------------------------------------------------------------#
# Add application targets
#------------------------------------------------------------------------------#
add_subdirectory("app/id_generators")
add_subdirectory("app/drivers")


#------------------------------------------------------------------------------#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#------------------------------------------------------------------------------#
