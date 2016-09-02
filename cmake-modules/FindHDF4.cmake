# Once done this will define
#  HDF4_FOUND - System has HDF4
#  HDF4_INCLUDE_DIRS - The HDF4 include directories
#  HDF4_LIBRARIES - The libraries needed to use HDF4
#  HDF4_DEFINITIONS - Compiler switches required for using HDF4

find_path(HDF4_INCLUDE_DIR hdf.h 
          HINTS ${PC_HDF4_INCLUDEDIR} ${PC_HDF4_INCLUDE_DIRS}
	  /usr/include/hdf
	  /usr/include
	  $ENV{HDF4_DIR}/include
          PATH_SUFFIXES)

find_library(HDF4_LIBRARY NAMES libmfhdf.a
             HINTS ${PC_HDF4_LIBDIR} 
		   ${PC_HDF4_LIBRARY_DIRS} 
	 	   /usr/lib 
		   $ENV{HDF4_DIR}/lib)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set HDF4_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(HDF4  DEFAULT_MSG
                                  HDF4_LIBRARY HDF4_INCLUDE_DIR)

mark_as_advanced(HDF4_INCLUDE_DIR HDF4_LIBRARY )

set(HDF4_LIBRARIES ${HDF4_LIBRARY} )
set(HDF4_INCLUDE_DIRS ${HDF4_INCLUDE_DIR} )
