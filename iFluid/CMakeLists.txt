file(GLOB ifluid_source "*.cpp" "*.h")
list(REMOVE_ITEM ifluid_source ${CMAKE_CURRENT_SOURCE_DIR}/iFluid.cpp)

add_library(ifluid STATIC ${ifluid_source})

target_link_libraries (ifluid frontier solver)

target_include_directories(ifluid PUBLIC 
				${frontier_include_dir}
                   		${solver_include_dir}
                              	${petsc_include_dir}
				${ifluid_include_dir}
                             	${external_include_dir})
