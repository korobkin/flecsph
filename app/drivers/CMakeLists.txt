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
# - SRC_PATH: location of main.cc, main_driver.cc for DRIVER_NAME
# - DIM_LIST: dimensions to compile DRIVER_NAME for, e.g. "1;2;3" for 1/2/3 dimension
function(add_driver driver_name src_path dim_list)
  foreach(dim ${dim_list})
    set(exe_name "${driver_name}_${dim}d")

    add_executable(${exe_name})

    target_sources(${exe_name}
      PRIVATE
        ${src_path}/main.cc
        ${src_path}/main_driver.cc
    )

    target_link_libraries(${exe_name}
      PRIVATE
        FleCSI::flecsi
        flecsph::library_flags
        flecsph::compile_flags
    )

    target_compile_definitions(${exe_name}
      PRIVATE
        EXT_GDIMENSION=${dim}
        FLECSI_ENABLE_SPECIALIZATION_TLT_INIT
        FLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER
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

  package_add_test(collapse_test test/collapse.cc newtonian/main_driver.cc)
  target_compile_options(collapse_test PRIVATE "-DEXT_GDIMENSION=3")
  configure_file(${PROJECT_SOURCE_DIR}/data/collapse_nx10.par ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

  add_custom_target(
      collapse_3d_generator_test
        COMMAND ${MPIEXEC} -n 1 ../id_generators/collapse_3d_generator ${PROJECT_SOURCE_DIR}/data/collapse_nx10.par
        DEPENDS ../id_generators/collapse_3d_generator
  )
  add_dependencies(collapse_test collapse_3d_generator_test)

  if(ENABLE_MPI_TESTS)

    package_add_test_MPI(collapse_MPI_test test/collapse.cc newtonian/main_driver.cc)
    target_compile_options(collapse_MPI_test PRIVATE "-DEXT_GDIMENSION=3")
    configure_file(${PROJECT_SOURCE_DIR}/data/collapse_nx10.par ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

    add_dependencies(collapse_MPI_test collapse_3d_generator_test)

  endif()  # ENABLE_MPI_TESTS
endif()  # ENABLE_UNIT_TESTS


# #------------------------------------------------------------------------------#
# # sodtube test, call the default parameter file
# #------------------------------------------------------------------------------#
if(ENABLE_UNIT_TESTS)

  package_add_test(sodtube_test test/sodtube.cc hydro/main_driver.cc)
  target_compile_options(sodtube_test PRIVATE "-DEXT_GDIMENSION=1")
  configure_file(${PROJECT_SOURCE_DIR}/data/sodtube_t1_n100.par ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

  add_custom_target(
    sodtube_1d_generator_test
      COMMAND ${MPIEXEC} -n 1 ../id_generators/sodtube_1d_generator ${PROJECT_SOURCE_DIR}/data/sodtube_t1_n100.par
      DEPENDS ../id_generators/sodtube_1d_generator
  )
  add_dependencies(sodtube_test sodtube_1d_generator_test)

  if(ENABLE_MPI_TESTS)

    package_add_test_MPI(sodtube_MPI_test test/sodtube.cc hydro/main_driver.cc)
    target_compile_options(sodtube_MPI_test PUBLIC "-DEXT_GDIMENSION=1")
    configure_file(${PROJECT_SOURCE_DIR}/data/sodtube_t1_n100.par ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
    add_dependencies(sodtube_MPI_test sodtube_1d_generator_test)
  endif()  # ENABLE_MPI_TESTS
endif()  # ENABLE_UNIT_TESTS
# #------------------------------------------------------------------------------#
# # RT test, call the default parameter file
# #------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)

  package_add_test(RT_test test/RT.cc hydro/main_driver.cc)
  target_compile_options(RT_test PRIVATE "-DEXT_GDIMENSION=2")
  configure_file(${PROJECT_SOURCE_DIR}/data/RT_2d.par ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

  add_custom_target(
    RT_2d_generator_test
      COMMAND ${MPIEXEC} -n 1 ../id_generators/RT_2d_generator ${PROJECT_SOURCE_DIR}/data/RT_2d.par
      DEPENDS ../id_generators/RT_2d_generator
  )
  add_dependencies(RT_test RT_2d_generator_test)
endif()

# #------------------------------------------------------------------------------#
# # KH test: uses KH_2d.par from data directory
# #------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)

  package_add_test(KH_test test/KH.cc hydro/main_driver.cc)
  target_compile_options(KH_test PRIVATE "-DEXT_GDIMENSION=2")
  configure_file(${PROJECT_SOURCE_DIR}/data/KH_2d.par ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

  add_custom_target(
    KH_2d_generator_test
      COMMAND ${MPIEXEC} -n 1 ../id_generators/KH_2d_generator ${PROJECT_SOURCE_DIR}/data/KH_2d.par
      DEPENDS ../id_generators/KH_2d_generator
  )
  add_dependencies(KH_test KH_2d_generator_test)
endif() # ENABLE_UNIT_TESTS

# #------------------------------------------------------------------------------#
# # sedov test with the default parameter file
# #------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)

  package_add_test(sedov_test test/sedov.cc hydro/main_driver.cc)
  target_compile_options(sedov_test PRIVATE "-DEXT_GDIMENSION=2")
  configure_file(${PROJECT_SOURCE_DIR}/data/sedov_nx20.par ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

  add_custom_target(
    sedov_2d_generator_test
      COMMAND ${MPIEXEC} -n 1 ../id_generators/sedov_2d_generator ${PROJECT_SOURCE_DIR}/data/sedov_nx20.par
      DEPENDS ../id_generators/sedov_2d_generator
  )
  add_dependencies(sedov_test sedov_2d_generator_test)
endif()

# #------------------------------------------------------------------------------#
# # relaxation test for the "mesa" potential in 3D using mesa_nx20.par file
# #------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)
  package_add_test(mesa_relaxation_test test/mesa.cc hydro/main_driver.cc)
  target_compile_options(mesa_relaxation_test PRIVATE "-DEXT_GDIMENSION=3")
  configure_file(${PROJECT_SOURCE_DIR}/data/mesa_nx20.par ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

  add_custom_target(
    mesa_3d_generator_test
      COMMAND ${MPIEXEC} -n 1 ../id_generators/sedov_3d_generator ${PROJECT_SOURCE_DIR}/data/mesa_nx20.par
      DEPENDS ../id_generators/sedov_3d_generator
  )
  add_dependencies(mesa_relaxation_test mesa_3d_generator_test)
endif()

# #------------------------------------------------------------------------------#
# # relaxation test for wvt in 3D using wvt_nx20.par file
# #------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)
  package_add_test(wvt_relaxation_test test/wvt.cc hydro/main_driver.cc)
  target_compile_options(wvt_relaxation_test PRIVATE "-DEXT_GDIMENSION=3")
  configure_file(${PROJECT_SOURCE_DIR}/data/wvt_nx20.par ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

  add_custom_target(
    wvt_3d_generator_test
      COMMAND ${MPIEXEC} -n 1 ../id_generators/sedov_3d_generator ${PROJECT_SOURCE_DIR}/data/wvt_nx20.par
      DEPENDS ../id_generators/sedov_3d_generator
  )
  add_dependencies(wvt_relaxation_test wvt_3d_generator_test)
endif()

# #------------------------------------------------------------------------------#
# # noh test, with the default parameter file
# #------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)
  package_add_test(noh_test test/noh.cc hydro/main_driver.cc)
  target_compile_options(noh_test PRIVATE "-DEXT_GDIMENSION=2")
  configure_file(${PROJECT_SOURCE_DIR}/data/noh_nx20.par ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

  add_custom_target(
    noh_2d_generator_test
      COMMAND ${MPIEXEC} -n 1 ../id_generators/noh_2d_generator ${PROJECT_SOURCE_DIR}/data/noh_nx20.par
      DEPENDS ../id_generators/noh_2d_generator
  )
  add_dependencies(noh_test noh_2d_generator_test)
endif()

# #------------------------------------------------------------------------------#
# # implosion test
# #------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)
  package_add_test(implosion_test test/implosion.cc hydro/main_driver.cc)
  target_compile_options(implosion_test PRIVATE "-DEXT_GDIMENSION=2")
  configure_file(${PROJECT_SOURCE_DIR}/data/implosion_nx20.par ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

  add_custom_target(
    implosion_2d_generator_test
      COMMAND ${MPIEXEC} -n 1 ../id_generators/implosion_2d_generator ${PROJECT_SOURCE_DIR}/data/implosion_nx20.par
      DEPENDS ../id_generators/implosion_2d_generator
  )
  add_dependencies(implosion_test implosion_2d_generator_test)
endif()
