file(GLOB solver_source "*.cpp" "*.h")

add_library(solver STATIC ${solver_source})

find_package (PETSc)

if (PETSc_FOUND)
    MESSAGE("-- Found PETSc: " ${PETSc_LIBRARY})
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__HYPRE__")
    list(APPEND external_include_dir ${PETSc_INCLUDE_DIRS})
else ()
    MESSAGE(FATAL_ERROR "Could NOT find PETSc, be sure to set PETSC_DIR and PETSc is configured with --download-hypre")
endif()


target_include_directories(solver PUBLIC 
                           	    ${frontier_include_dir} 
				    ${solver_inlcude_dir} 
				    ${petsc_include_dir}
				    ${external_include_dir})
target_link_libraries(solver frontier ${PETSc_LIBRARY})
