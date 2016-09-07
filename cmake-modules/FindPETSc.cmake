# Once done this will define
#  PETSc_FOUND - System has PETSc
#  PETSc_INCLUDE_DIRS - The PETSc include directories
#  PETSc_LIBRARIES - The libraries needed to use PETSc
#  PETSc_DEFINITIONS - Compiler switches required for using PETSc

find_path(PETSc_INCLUDE_DIR petsc.h 
          HINTS
	  /usr/include
	  /usr/local/include
	  $ENV{PETSC_DIR}/include
          PATH_SUFFIXES petsc/include)

find_library(PETSc_LIBRARY NAMES petsc
             HINTS $ENV{PETSC_DIR}/lib)

find_library(HYPRE_LIBRARY NAMES HYPRE
             HINTS $ENV{PETSC_DIR}/lib)

set (PETSc_LIBRARY ${PETSc_LIBRARY} ${HYPRE_LIBRARY})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PETSc_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(PETSc  DEFAULT_MSG
                                  PETSc_LIBRARY PETSc_INCLUDE_DIR)

mark_as_advanced(PETSc_INCLUDE_DIR PETSc_LIBRARY )

set(PETSc_LIBRARIES ${PETSc_LIBRARY} )
set(PETSc_INCLUDE_DIRS ${PETSc_INCLUDE_DIR} )
