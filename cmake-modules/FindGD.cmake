# Once done this will define
#  GD_FOUND - System has GD
#  GD_INCLUDE_DIRS - The GD include directories
#  GD_LIBRARIES - The libraries needed to use GD
#  GD_DEFINITIONS - Compiler switches required for using GD

find_path(GD_INCLUDE_DIR gd.h 
          HINTS ${PC_GD_INCLUDEDIR} ${PC_GD_INCLUDE_DIRS}
	  $ENV{GD_DIR}/include
	  /usr/include
	  /usr/local/include
          PATH_SUFFIXES)

find_library(GD_LIBRARY NAMES libgd.a 
             HINTS ${PC_GD_LIBDIR} ${PC_GD_LIBRARY_DIRS} 
		$ENV{GD_DIR}/lib
		/usr/lib)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GD_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GD  DEFAULT_MSG
                                  GD_LIBRARY GD_INCLUDE_DIR)

mark_as_advanced(GD_INCLUDE_DIR GD_LIBRARY )

set(GD_LIBRARIES ${GD_LIBRARY} )
set(GD_INCLUDE_DIRS ${GD_INCLUDE_DIR} )
