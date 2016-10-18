cmake_minimum_required (VERSION 2.8.11)
project (climate C CXX)

set (CMAKE_C_FLAGS)
set (CMAKE_CXX_FLAGS "-fopenmp")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -rdynamic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -rdynamic")

#detect dependencies
list (APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")

find_package (MPI REQUIRED)

find_package (HDF4 QUIET)

find_package (GD QUIET)

if (MPI_FOUND) 
    MESSAGE(STATUS "Found MPI: ${MPI_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__MPI__")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__MPI__")
else ()
    MESSAGE(FATAL_ERROR "Could NOT find MPI")
endif()

if (HDF4_FOUND)
    MESSAGE(STATUS "Found HDF4: ${HDF4_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_HDF")
else ()
    MESSAGE(WARNING "HDF4 cannot be found. Set HDF4_ROOT is hdf animation is needed. Build continues ...")
endif()

if (GD_FOUND)
    MESSAGE(STATUS "Found GD: ${GD_LIBRARY}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__GD__")
else ()
    MESSAGE(WARNING "GD cannot be found. Set GD_DIR if gd animation is needed. Build continues ...")
endif()

#dependencies detection finished

set(frontier_include_dir ${CMAKE_SOURCE_DIR}/src/front
                         ${CMAKE_SOURCE_DIR}/src/intfc
                         ${CMAKE_SOURCE_DIR}/src/util
                         ${CMAKE_SOURCE_DIR}/src)

set(solver_include_dir ${CMAKE_SOURCE_DIR}/solver)

set(ifluid_include_dir ${CMAKE_SOURCE_DIR}/iFluid)

set(climate_include_dir ${CMAKE_SOURCE_DIR}/climate)

set(project_include_dir ${frontier_include_dir}
			${solver_include_dir}
			${ifluid_include_dir}
			${climate_include_dir})

#external libraries include directories
set(external_include_dir 
		${MPI_C_INCLUDE_PATH})

#external libraries link dir
set(external_lib_dir 
                ${MPI_C_LIBRARIES}) 

if (GD_FOUND)
	list( APPEND external_include_dir ${GD_INCLUDE_DIR})
	list( APPEND external_lib_dir ${GD_LIBRARY})
endif()

if (HDF4_FOUND)
	list( APPEND external_include_dir ${HDF4_INCLUDE_DIR})
	list( APPEND external_lib_dir ${HDF4_LIBRARIES})
endif()

# Recurse into the following subdirectories
add_subdirectory (climate)
add_subdirectory (iFluid)
add_subdirectory (src)
add_subdirectory (solver)
