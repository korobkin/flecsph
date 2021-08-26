#------------------------------------------------------------------------------#
# Find FleCSI, This should be in FindFleCSI
#------------------------------------------------------------------------------#
FIND_PATH(FleCSI_INCLUDE_DIR NAMES flecsi.h flecsi-config.h)
# Look for the library.
FIND_LIBRARY(FleCSI_LIBRARY NAMES flecsi libflecsi FleCSI libFleCSI)

# handle the QUIETLY and REQUIRED arguments and set FleCSI_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FleCSI DEFAULT_MSG FleCSI_INCLUDE_DIR FleCSI_LIBRARY)

if(FleCSI_FOUND)
    # Copy the results to the output variables.
    SET(FleCSI_INCLUDE_DIRS ${FleCSI_INCLUDE_DIR})
    SET(FleCSI_LIBRARIES ${FleCSI_LIBRARY})
else()
    SET(FleCSI_INCLUDE_DIRS)
    SET(FleCSI_LIBRARIES)
endif()

add_library(FleCSI::flecsi UNKNOWN IMPORTED)
set_target_properties(FleCSI::flecsi PROPERTIES
    IMPORTED_LOCATION "${FleCSI_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FleCSI_INCLUDE_DIR}"
)

MARK_AS_ADVANCED(FleCSI_INCLUDE_DIR FleCSI_LIBRARY)

