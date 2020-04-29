#------------------------------------------------------------------------------#
# Copyright (c) 2020 Triad National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

project(flecsph)
# TODO: is this still necessary?
cmake_minimum_required(VERSION 3.12)

#------------------------------------------------------------------------------#
# The project name
#------------------------------------------------------------------------------#
project(flecsph)

option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
option(ENABLE_DEBUG_TREE "Enable debug tree" OFF)
option(ENABLE_DEBUG "Compile in DEBUG mode" OFF)
set(LOG_STRIP_LEVEL 0 CACHE STRING "LOG Strip level")
set(KEY_INTEGER_TYPE "uint64_t" CACHE STRING "Type of integer used to generate keys")
set_property(CACHE KEY_INTEGER_TYPE PROPERTY STRINGS "uint32_t" "uint64_t" "uint128_t")

# Set a default build type
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# create directory if tests 
if(ENABLE_UNIT_TESTS)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
endif()

# interface libraries
# see: cmake/SetupDeps.cmake, cmake/SetupFlags.cmake
add_library(flecsph::compile_flags INTERFACE IMPORTED)
add_library(flecsph::library_flags INTERFACE IMPORTED)

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(SetupDeps)
include(SetupFlags)
include(SetupTesting)

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
