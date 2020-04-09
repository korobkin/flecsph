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
# this may, in future, be better placed in it's own module
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(flecsph::compile_flags
        INTERFACE
            "$<$<CONFIG:Debug>:-g>"
            "$<$<CONFIG:Release>:-Ofast;-march=native;-mtune=native>"
    )
    # target_link_options(flecsph::compile_flags
    #     INTERFACE
    #         "$<$<CONFIG:Release>:-flto>"
    # )
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    target_compile_options(flecsph::compile_flags
        INTERFACE
            "$<$<CONFIG:Debug>:-g;-traceback>"
            "$<$<CONFIG:Release>:-fast;-xHost>"
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
# NOTE: imported libraries bring in includes, definitions, libs, &t. convienent!
target_link_libraries(flecsph::library_flags
    INTERFACE
        Threads::Threads
        OpenMP::OpenMP_CXX
        MPI::MPI_CXX
        GSL::gsl
        Boost::headers
        m
        "$<$<BOOL:${ENABLE_UNIT_TESTS}>:GTest::GTest>"
        "$<$<BOOL:${ENABLE_UNIT_TESTS}>:GTest::Main>"
        ${HDF5_LIBRARIES}
)

# HDF5 doesn't provide imported interface (as far as I can tell),
# so explicitily provide
target_include_directories(flecsph::compile_flags
    INTERFACE
        ${HDF5_INCLUDE_DIR}
)