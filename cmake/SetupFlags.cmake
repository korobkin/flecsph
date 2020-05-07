# set C++17
target_compile_features(flecsph::compile_flags
    INTERFACE
        cxx_std_17)

# some general definitions; subdirectories may define
# targets with specialized definitions
target_compile_definitions(flecsph::compile_flags
    INTERFACE
        "LOG_STRIP_LEVEL=${LOG_STRIP_LEVEL}"
        "$<$<BOOL:${ENABLE_DEBUG_TREE}>:ENABLE_DEBUG_TREE=${ENABLE_DEBUG_TREE}>"
        "PARALLEL_IO"
)

# compiler-specific flags
# TODO: future, may be moved to it's own module
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(flecsph::compile_flags
        INTERFACE
            "$<$<CONFIG:Debug>:-O2>"
            "$<$<CONFIG:Release>:-floop-nest-optimize>"
    )
    # TODO: Check if LTO is worth doing
    # target_link_options(flecsph::compile_flags
    #     INTERFACE
    #         "$<$<CONFIG:Release>:-flto>"
    # )
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Cray")
    target_compile_options(flecsph::compiler_flags
      INTERFACE
        "$<$<CONFIG:Debug:-O2>"
    )
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    target_compile_options(flecsph::compile_flags
        INTERFACE
            "$<$<CONFIG:Debug>:-g;-O2;-traceback>"
            "$<$<CONFIG:Release>:-O3>"
    )
    target_link_options(flecsph::library_flags
        INTERFACE
            # mpi issues with ipo
            "-no-ipo"
    )
endif()

# global includes
target_include_directories(flecsph::compile_flags
    INTERFACE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/physics
        ${CMAKE_SOURCE_DIR}/include/physics/eos
        ${CMAKE_SOURCE_DIR}/include/physics/gw_rad
        ${CMAKE_SOURCE_DIR}/app/drivers/include
        ${CMAKE_SOURCE_DIR}/mpisph
)

# global libraries
# NOTE: imported libraries bring in includes, definitions, libs; convienent!
target_link_libraries(flecsph::library_flags
    INTERFACE
        Threads::Threads
        OpenMP::OpenMP_CXX
        MPI::MPI_CXX
        GSL::gsl
        Boost::headers
        m
        "$<$<BOOL:${ENABLE_UNIT_TESTS}>:"
          "GTest::GTest"
          "GTest::Main"
        ">"
        ${HDF5_LIBRARIES}
)

# HDF5 doesn't provide imported interface (as far as I can tell),
# so explicitily provide
target_include_directories(flecsph::compile_flags
    INTERFACE
        ${HDF5_INCLUDE_DIR}
)
