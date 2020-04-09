#----------------------------------------------------------------------------#
# Copyright (c) 2017 Triad National Security, LLC
# All rights reserved.
#----------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Executables
#------------------------------------------------------------------------------#

# see 'drivers/CMakeLists.txt` for explaination
# very similar function to driver, may consider
# generalizing
function(add_generator generator_name dim_list)
  foreach(dim ${dim_list})
    set(exe_name "${generator_name}_${dim}d_generator")

    add_executable(${exe_name})

    target_sources(${exe_name}
      PRIVATE
        ${generator_name}/main.cc
    )

    target_compile_definitions(${exe_name}
      PRIVATE
        EXT_GDIMENSION=${dim}
    )

    # NOTE: using imported 'FleCSI::flecsi' will
    # also append runtime_driver.cc, causing unresolved references
    # at linktime.
    # Use the explicit include to avoid this
    target_include_directories(${exe_name}
      PRIVATE
        ${CMAKE_SOURCE_DIR}/app/id_generators/include
        ${FleCSI_INCLUDE_DIR}
    )

    target_link_libraries(${exe_name}
      PRIVATE
        flecsph::library_flags
        flecsph::compile_flags
    )

    install(TARGETS ${exe_name} RUNTIME DESTINATION bin)
  endforeach()
endfunction()

add_generator(sodtube "1;2;3")
add_generator(sedov "1;2;3")
add_generator(noh "1;2;3")
add_generator(implosion "2;3")
add_generator(collapse "3")
add_generator(wtunnel "2;3")
add_generator(KH "2;3")
add_generator(RT "2;3")
add_generator(fluid "2;3")
