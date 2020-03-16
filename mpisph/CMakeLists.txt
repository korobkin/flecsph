#~----------------------------------------------------------------------------~#
#  Copyright (c) 2017 Triad National Security, LLC
#  All rights reserved.
#~----------------------------------------------------------------------------~#

#~----------------------------------------------------------------------------~#
#
# /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
# /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
# /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
# /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
# /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
# /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
# /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
# //       ///  //////   //////  ////////  //       //      //
#
#~---------------------------------------------------------------------------~*#

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/bodies/)
include_directories(${CMAKE_SOURCE_DIR}/include/physics/default/)
include_directories(${CMAKE_SOURCE_DIR}/include/physics/)
include_directories(${CMAKE_SOURCE_DIR}/include/physics/eos)
include_directories(${CMAKE_SOURCE_DIR}/mpisph/test)
include_directories(${CMAKE_SOURCE_DIR}/mpisph/)


#------------------------------------------------------------------------------#
# tree topology tests
#------------------------------------------------------------------------------#

package_add_test(tree3d test/tree3d.cc)
package_add_test(mpi_qsort test/mpi_qsort.cc)

package_add_test(io test/io.cc)
configure_file(test/io_test.h5part ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

package_add_test(bs test/bs.cc)
configure_file(test/io_test.h5part ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

#~---------------------------------------------------------------------------~-#
# Formatting options
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
