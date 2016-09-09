# Once done this will define
#  GD_FOUND - System has GD
#  GD_INCLUDE_DIRS - The GD include directories
#  GD_LIBRARIES - The libraries needed to use GD
#  GD_DEFINITIONS - Compiler switches required for using GD
find_program(GD2TOGIF gd2togif HINTS $ENV{PATH})

if (NOT GD2TOGIF)
    set(GD_DIR, "")
else()
    get_filename_component(GD_DIR ${GD2TOGIF}/../.. ABSOLUTE)
endif()

find_path(GD_INCLUDE_DIR gd.h 
          HINTS ${GD_DIR}/include
		/usr/include
		/usr/local/include)

find_library(GD_LIBRARY NAMES libgd.a 
             HINTS ${GD_DIR}/lib
		   /usr/lib
		   /usr/lib64
		   /usr/local/lib)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GD_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GD  DEFAULT_MSG
                                  GD_LIBRARY GD_INCLUDE_DIR)

mark_as_advanced(GD_INCLUDE_DIR GD_LIBRARY )

set(GD_LIBRARIES ${GD_LIBRARY} )
set(GD_INCLUDE_DIRS ${GD_INCLUDE_DIR} )
