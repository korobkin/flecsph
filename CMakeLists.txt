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

set(ENABLE_UNIT_TESTS ON CACHE BOOL "Enable unit tests")
set(ENABLE_DEBUG_TREE OFF CACHE BOOL "Enable debug tree")
set(ENABLE_DEBUG OFF CACHE BOOL "Compile in DEBUG mode" )
set(LOG_STRIP_LEVEL 0 CACHE STRING "LOG Strip level")

if(ENABLE_DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLOG_STRIP_LEVEL=${LOG_STRIP_LEVEL}")

if(ENABLE_DEBUG_TREE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_DEBUG_TREE=${ENABLE_DEBUG_TREE}")
endif()

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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17 -Wno-return-type")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wno-return-type")
endif()

#------------------------------------------------------------------------------#
# Find FleCSI, This should be in FindFleCSI
#------------------------------------------------------------------------------#
FIND_PATH(FleCSI_INCLUDE_DIR NAMES flecsi.h flecsi-config.h)
# Look for the library.
FIND_LIBRARY(FleCSI_LIBRARY NAMES flecsi libflecsi FleCSI libFleCSI)

# Look for the runtime driver
FIND_PATH(FleCSI_RUNTIME
	NAMES
		runtime_driver.cc
	PATH_SUFFIXES
		share/FleCSI/runtime share/flecsi/runtime
)

# handle the QUIETLY and REQUIRED arguments and set FleCSI_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FleCSI DEFAULT_MSG FleCSI_INCLUDE_DIR FleCSI_LIBRARY FleCSI_RUNTIME)

# Copy the results to the output variables.
SET(FleCSI_INCLUDE_DIRS ${FleCSI_INCLUDE_DIR})
SET(FleCSI_LIBRARIES ${FleCSI_LIBRARY})
SET(FleCSI_RUNTIME ${FleCSI_RUNTIME})

MARK_AS_ADVANCED(FleCSI_INCLUDE_DIR FleCSI_LIBRARY FleCSI_RUNTIME)

set(FleCSPH_LIBRARIES ${FleCSI_LIBRARIES})

#------------------------------------------------------------------------------#
# Find GoogleTest
#------------------------------------------------------------------------------#
if(ENABLE_UNIT_TESTS)
enable_testing()
find_package(GTest REQUIRED)
list (APPEND FleCSPH_LIBRARIES GTest::GTest GTest::Main)
endif()
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
