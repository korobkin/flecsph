#----------------------------------------------------------------------------#
# Copyright (c) 2017 Triad National Security, LLC
# All rights reserved.
#----------------------------------------------------------------------------#
option(ENABLE_MPI_TESTS "Enable unit testing with MPI" ON)

#------------------------------------------------------------------------------#
# Debug and release flags
#------------------------------------------------------------------------------#

# UPDATE mauneyc, 04-04-2020
# flags set in "cmake/SetupFlags.cmake"
# keeping these here for reference if needed
# For vectorization: -fopt-info-vec-missed/optimized/note/all
# set(CMAKE_CXX_FLAGS_DEBUG
#   "-mavx2 -msse4 -ftree-vectorize -march=native \
#   -pthread -lpthread -DPARALLEL_IO -Wno-sign-compare -Wno-reorder -Wno-narrowing \
#   -Wno-deprecated-declarations  -ffast-math -g -O3 -Wall -Wno-return-type \
# -Wno-unused -Wno-comment -Wno-parentheses -lm")
# set(CMAKE_CXX_FLAGS_RELEASE
#  "-march=native -pthread -lpthread -DPARALLEL_IO -Wno-sign-compare -Wno-reorder \
#   -Wno-narrowing -Wno-deprecated-declarations -ftree-vectorize -ffast-math -msse2 \
#   -O3 -Wall -Wno-return-type -Wno-unused -Wno-comment -Wno-parentheses")

# adds a driver executble
# - DRIVER_NAME: executable name of the driver
# - SRC_PATH: location of main_driver.cc for DRIVER_NAME
# - DIM_LIST: dimensions to compile DRIVER_NAME for, e.g. "1;2;3" for 1/2/3 dimension
function(add_driver driver_name src_path dim_list)
  foreach(dim ${dim_list})
    set(exe_name "${driver_name}_${dim}d")

    add_executable(${exe_name})

    target_sources(${exe_name}
      PRIVATE
        ${src_path}/main_driver.cc
    )

    target_link_libraries(${exe_name}
      PRIVATE
        flecsph::flags
        FleCSI::flecsi
    )

    target_compile_definitions(${exe_name}
      PRIVATE
        "EXT_GDIMENSION=${dim}"
        "FLECSI_ENABLE_SPECIALIZATION_TLT_INIT"
        "FLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER"
    )
    install(TARGETS ${exe_name} RUNTIME DESTINATION bin)

  endforeach()
endfunction()
#------------------------------------------------------------------------------#
# Hydro drivers without gravity
#------------------------------------------------------------------------------#
add_driver(hydro hydro "1;2;3")

# #------------------------------------------------------------------------------#
# # Tree drivers
# #------------------------------------------------------------------------------#
add_driver(tree tree "1;2;3")

# #------------------------------------------------------------------------------#
# # WVT drivers
# #------------------------------------------------------------------------------#
add_driver(wvt wvt "2;3")

# #------------------------------------------------------------------------------#
# # Hydro drivers with Newtonian gravity
# #------------------------------------------------------------------------------#
add_driver(newtonian newtonian "3")

# #------------------------------------------------------------------------------#
# # collapse test, call the default parameter file
# #------------------------------------------------------------------------------#

# these could be function'd, but not sure
if(ENABLE_UNIT_TESTS)

  package_add_test(collapse_test collapse_nx10.par newtonian_3d)
  set_tests_properties(
    collapse_test PROPERTIES DEPENDS collapse_3d_generator_collapse_test)
  set_tests_properties(collapse_test PROPERTIES DEPENDS newtonian_3d)

  if(ENABLE_MPI_TESTS)

    package_add_test_MPI(collapse_MPI_test collapse_nx10.par newtonian_3d)
    set_tests_properties(
      collapse_MPI_test PROPERTIES DEPENDS collapse_3d_generator_collapse_test)
    set_tests_properties( collapse_MPI_test PROPERTIES DEPENDS newtonian_3d)

  endif()  # ENABLE_MPI_TESTS

# #------------------------------------------------------------------------------#
# # sodtube test, call the default parameter file
# #------------------------------------------------------------------------------#

  package_add_test(sodtube_test sodtube_t1_n100.par hydro_1d)
  set_tests_properties( sodtube_test PROPERTIES DEPENDS sodtube_1d_generator_sodtube_test)
  set_tests_properties( sodtube_test PROPERTIES DEPENDS hydro_1d)

  if(ENABLE_MPI_TESTS)

    package_add_test_MPI(sodtube_MPI_test sodtube_t1_n100.par hydro_1d)
    set_tests_properties(sodtube_MPI_test PROPERTIES DEPENDS sodtube_1d_generator_sodtube_test)
    set_tests_properties(sodtube_MPI_test PROPERTIES DEPENDS hydro_1d)

  endif()  # ENABLE_MPI_TESTS
# #------------------------------------------------------------------------------#
# # RT test, call the default parameter file
# #------------------------------------------------------------------------------#

  package_add_test(RT_test RT_2d.par hydro_2d)
  set_tests_properties(RT_test PROPERTIES DEPENDS RT_2d_generator_RT_test)
  set_tests_properties(RT_test PROPERTIES DEPENDS hydro_2d)

# #------------------------------------------------------------------------------#
# # KH test: uses KH_2d.par from data directory
# #------------------------------------------------------------------------------#

  package_add_test(KH_test KH_2d.par hydro_2d)
  set_tests_properties( KH_test PROPERTIES DEPENDS KH_2d_generator_KH_test)
  set_tests_properties( KH_test PROPERTIES DEPENDS hydro_2d)

# #------------------------------------------------------------------------------#
# # sedov test with the default parameter file
# #------------------------------------------------------------------------------#

  package_add_test(sedov_test sedov_nx20.par hydro_2d)
  set_tests_properties( sedov_test PROPERTIES DEPENDS sedov_2d_generator_sedov_test)
  set_tests_properties( sedov_test PROPERTIES DEPENDS hydro_2d)

# #------------------------------------------------------------------------------#
# # relaxation test for the "mesa" potential in 3D using mesa_nx20.par file
# #------------------------------------------------------------------------------#
  package_add_test(mesa_relaxation_test mesa_nx20.par hydro_3d)
  set_tests_properties( mesa_relaxation_test PROPERTIES DEPENDS sedov_3d_generator_mesa_test)
  set_tests_properties( mesa_relaxation_test PROPERTIES DEPENDS hydro_3d)

# #------------------------------------------------------------------------------#
# # relaxation test for wvt in 3D using wvt_nx20.par file
# #------------------------------------------------------------------------------#

  package_add_test(wvt_relaxation_test wvt_nx20.par hydro_3d)
  set_tests_properties( wvt_relaxation_test PROPERTIES DEPENDS sedov_3d_generator_wvt_test)
  set_tests_properties( wvt_relaxation_test PROPERTIES DEPENDS hydro_3d)


# #------------------------------------------------------------------------------#
# # noh test, with the default parameter file
# #------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)
  package_add_test(noh_test noh_nx20.par hydro_2d)
  set_tests_properties( noh_test PROPERTIES DEPENDS noh_2d_generator_noh_test)
  set_tests_properties( noh_test PROPERTIES DEPENDS hydro_2d)

# #------------------------------------------------------------------------------#
# # implosion test
# #------------------------------------------------------------------------------#

  package_add_test(implosion_test implosion_nx20.par hydro_2d)
  set_tests_properties( implosion_test PROPERTIES DEPENDS implosion_2d_generator_implositon_test)
  set_tests_properties( implosion_test PROPERTIES DEPENDS hydro_2d)


endif()
endif()
