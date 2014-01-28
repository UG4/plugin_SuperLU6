################################################################################
# SuperLU
################################################################################

cmake_minimum_required(VERSION 2.6)

project(SuperLU)

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")


# sources
set(srcSuperLU 	super_lu_bridge.cpp super_lu.cpp super_lu.h)

# if super lu is ON
if(SUPERLU_PATH)   
    set(INTERNAL_SUPERLU_PATH "${SUPERLU_PATH}")
    message(STATUS "Info: SuperLU is ON, using user provided path ${INTERNAL_SUPERLU_PATH}")
else(SUPERLU_PATH)
    set(INTERNAL_SUPERLU_PATH "$ENV{HOME}/local/SuperLU/used/")
    message(STATUS "Info: SuperLU is ON, using uginstall provided path ${INTERNAL_SUPERLU_PATH}")    
endif(SUPERLU_PATH)


include_directories("${INTERNAL_SUPERLU_PATH}/SRC/")


add_library(SuperLULib STATIC IMPORTED)
set_property(TARGET SuperLULib PROPERTY IMPORTED_LOCATION "${INTERNAL_SUPERLU_PATH}/lib/libsuperlu_4.3.a")
 
if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/superlu ${srcSuperLU})
	target_link_libraries(ug4 SuperLULib)
else(buildEmbeddedPlugins)
	add_library(superLU SHARED ${srcSuperLU})	
	target_link_libraries (superLU ug4 )
	target_link_libraries(superLU SuperLULib)
endif(buildEmbeddedPlugins)
