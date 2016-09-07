file(GLOB climate_source "*.cpp" "*.h")

#create excutable
add_executable(climate ${climate_source})

#this project needs fftw3 package
find_package (FFTW3)

if (FFTW3_FOUND)
    MESSAGE("-- Found FFTW3: " ${FFTW3_LIBRARY})
else ()
    MESSAGE(FATAL_ERROR "FFTW3 cannot be found. Be sure to set FFTW3_DIR.")
endif()

#find hdf5 library for big data
find_path(PETSC_HDF5_PATH petscviewerhdf5.h 
		PATHS $ENV{PETSC_DIR}/include)

if (NOT ${PETSC_HDF5_PATH} STREQUAL "")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__HDF5__")
    MESSAGE(STATUS "PETSc was compiled with HDF5. hdf5 data format is enabled")
else ()
    MESSAGE(WARNING "PETSc is not configured with hdf5. Configure PETSc with --download-hdf5 to restore raw data")
endif()

#end find hdf5 library

if (FFTW3_FOUND)
        list( APPEND external_include_dir ${FFTW3_INCLUDE_DIR})
        list( APPEND external_lib_dir ${FFTW3_LIBRARY})
endif()

#include directories
target_include_directories(climate PUBLIC
				   ${solver_include_dir}
				   ${frontier_include_dir}	
				   ${ifluid_include_dir}
				   ${external_include_dir})

#link libraries
target_link_libraries(climate frontier solver ifluid ${FFTW3_LIBRARY} ${MPI_C_LIBRARIES})
