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
