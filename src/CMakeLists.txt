file(GLOB frontier_source 
	"front/*.h"
	"front/*.c"
	"intfc/*.h"
	"intfc/*.c"
	"util/*.h"
	"util/*.c"
	"FronTier.h"
)
list(REMOVE_ITEM frontier_source ${CMAKE_CURRENT_SOURCE_DIR}/intfc/testintfc.c)
set(frontier_dir ${CMAKE_CURRENT_SOURCE_DIR})

add_library(frontier STATIC ${frontier_source}) 
target_include_directories(frontier PUBLIC 
				${frontier_include_dir} 
				${external_include_dir})
set(link_libs)
if (HDF4_FOUND)
    list (APPEND link_libs ${HDF4_LIBRARIES})
endif()

if (GD_FOUND)
    list (APPEND link_libs ${GD_LIBRARY})
endif()

target_link_libraries(frontier ${link_libs})
