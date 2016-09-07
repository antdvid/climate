#
# To be used by projects that make use of CMakeified hdf-4.2.12
#

#
# Find the HDF4 includes and get all installed hdf4 library settings from
# HDF4-config.cmake file : Requires a CMake compatible hdf-4.2.12 or later 
# for this feature to work. The following vars are set if hdf4 is found.
#
# HDF4_FOUND               - True if found, otherwise all other vars are undefined
# HDF4_INCLUDE_DIR         - The include dir for main *.h files
# HDF4_FORTRAN_INCLUDE_DIR - The include dir for fortran modules and headers
# HDF4_VERSION_STRING      - full version (e.g. 4.2.12)
# HDF4_VERSION_MAJOR       - major part of version (e.g. 4.2)
# HDF4_VERSION_MINOR       - minor part (e.g. 12)
# 
# The following boolean vars will be defined
# HDF4_ENABLE_PARALLEL - 1 if HDF4 parallel supported
# HDF4_BUILD_FORTRAN   - 1 if HDF4 was compiled with fortran on
# HDF4_BUILD_CPP_LIB   - 1 if HDF4 was compiled with cpp on
# HDF4_BUILD_TOOLS     - 1 if HDF4 was compiled with tools on
# 
# Target names that are valid (depending on enabled options)
# will be the following
#
# hdf              : HDF4 C library
# hdf_f90cstub     : used by Fortran to C interface
# hdf_fortran      : Fortran HDF4 library
# mfhdf            : HDF4 multi-file C interface library
# xdr              : RPC library
# mfhdf_f90cstub   : used by Fortran to C interface to multi-file library
# mfhdf_fortran    : Fortran multi-file library
# 
# To aid in finding HDF4 as part of a subproject set
# HDF4_ROOT_DIR_HINT to the location where hdf4-config.cmake lies

include (SelectLibraryConfigurations)
include (FindPackageHandleStandardArgs)

# The HINTS option should only be used for values computed from the system.
set (_HDF4_HINTS
    $ENV{HOME}/.local
    $ENV{HDF4_ROOT}
    $ENV{HDF4_ROOT_DIR_HINT}
)
# Hard-coded guesses should still go in PATHS. This ensures that the user
# environment can always override hard guesses.
set (_HDF4_PATHS
    $ENV{HOME}/.local
    $ENV{HDF4_ROOT}
    $ENV{HDF4_ROOT_DIR_HINT}
    /usr/lib/hdf4
    /usr/share/hdf4
    /usr/local/hdf4
    /usr/local/hdf4/share
)

FIND_PATH (HDF4_ROOT_DIR "hdf4-config.cmake"
    HINTS ${_HDF4_HINTS}
    PATHS ${_HDF4_PATHS}
    PATH_SUFFIXES
        cmake/hdf4
        lib/cmake/hdf4
        share/cmake/hdf4
)

FIND_PATH (HDF4_INCLUDE_DIRS "hdf.h"
    HINTS ${_HDF4_HINTS}
    PATHS ${_HDF4_PATHS}
    PATH_SUFFIXES
        include
        Include
)

# Find hdf4 library and its dependencies, modified by zheng.gao@stonybrook.edu
# Search path is under HDF4_ROOT and the same level in case other libraries
# are installed with HDF4
if ( "$ENV{HDF4_ROOT}" STREQUAL "" )
    MESSAGE(WARNING "HDF4_ROOT is not set, HDF4 is not used")
    set(HDF4_PARENT_DIR)
else()
    get_filename_component(HDF4_PARENT_DIR $ENV{HDF4_ROOT} DIRECTORY)
endif()
set(HDF_SEARCH_PATH /usr/lib
		    $ENV{HDF4_ROOT}/lib
		    ${HDF4_PARENT_DIR})

FIND_LIBRARY(HDF4_LIBRARIES NAMES mfhdf HINTS ${HDF_SEARCH_PATH} PATH_SUFFIXES hdf/lib hdf4/lib)
FIND_LIBRARY(JPEG_LIBRARIES NAMES jpeg  HINTS ${HDF_SEARCH_PATH} PATH_SUFFIXES jpeg/lib)
FIND_LIBRARY(DF_LIBRARIES   NAMES df    HINTS ${HDF_SEARCH_PATH} PATH_SUFFIXES df/lib) 
FIND_LIBRARY(SZIP_LIBRARIES NAMES sz  HINTS ${HDF_SEARCH_PATH} PATH_SUFFIXES szip/lib)
FIND_LIBRARY(Z_LIBRARIES    NAMES z     HINTS ${HDF_SEARCH_PATH} PATH_SUFFIXES zlib/lib z/lib)
LIST(APPEND HDF4_LIBRARIES ${JPEG_LIBRARIES}
			   ${DF_LIBRARIES}
			   ${SZIP_LIBRARIES}
			   ${Z_LIBRARIES})

# For backwards compatibility we set HDF4_INCLUDE_DIR to the value of
# HDF4_INCLUDE_DIRS
set ( HDF4_INCLUDE_DIR "${HDF4_INCLUDE_DIRS}" )

if (HDF4_INCLUDE_DIR)
  set (HDF4_FOUND "YES")
endif (HDF4_INCLUDE_DIR)
