################################################################################
#      SuperLU plugin
#----------------------------
# information about SuperLU is here http://crd-legacy.lbl.gov/~xiaoye/SuperLU
# if you are using SuperLU, remember to CITE it!
# bibtex information is here
#     http://crd-legacy.lbl.gov/~xiaoye/SuperLU/referencing.html
#
# ug4 plugin by Martin Rupp
# how to use:
#     cmake -DSuperLU=ON ..
# this will install SuperLU with uginstall to ~/local/SuperLU
# this is tested for Macs and Cekon
# if this does not work, you can supply your own SuperLU path with
#     cmake -DSuperLU=ON -DSUPERLU_PATH=/path/to/my/superlu .. 
################################################################################

cmake_minimum_required(VERSION 2.6)

project(SuperLU)

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")


# sources
set(srcSuperLU 	super_lu_bridge.cpp super_lu.cpp super_lu.h)

# detect path for SuperLU-Library
if(SUPERLU_PATH)
    # the user can suply a SUPER_LU path for his own compiled SuperLU   
    set(INTERNAL_SUPERLU_PATH "${SUPERLU_PATH}")
    message(STATUS "Info: SuperLU is ON, using user provided path ${INTERNAL_SUPERLU_PATH}")
else(SUPERLU_PATH)
    # in other cases we use uginstall to download & install superlu (works on machines like Mac and Cekon)
    if(NOT EXISTS "$ENV{HOME}/local/SuperLU/used/installed")
        execute_process(COMMAND bash ${UG_ROOT_PATH}/scripts/shell/uginstall SuperLU)
        if(NOT EXISTS "$ENV{HOME}/local/SuperLU/used/installed")
            message(ERROR "SuperLU could not be installed by uginstall.")
        endif(NOT EXISTS "$ENV{HOME}/local/SuperLU/used/installed")        
    endif(NOT EXISTS "$ENV{HOME}/local/SuperLU/used/installed")
    
    set(INTERNAL_SUPERLU_PATH "$ENV{HOME}/local/SuperLU/used/")
    
    message(STATUS "Info: SuperLU is ON, using uginstall provided path ${INTERNAL_SUPERLU_PATH}")    
    
endif(SUPERLU_PATH)

# Add to build process 
if(buildEmbeddedPlugins)
	ExportSources(plugins/experimental/superlu ${srcSuperLU})
    ExportDependencies("${INTERNAL_SUPERLU_PATH}/lib/libsuperlu_4.3.a")
    if(EXISTS "${INTERNAL_SUPERLU_PATH}lib/blas.a")
        ExportDependencies("${INTERNAL_SUPERLU_PATH}/lib/blas.a")
    endif(EXISTS "${INTERNAL_SUPERLU_PATH}lib/blas.a")
    ExportIncludes("${INTERNAL_SUPERLU_PATH}/SRC/")
else(buildEmbeddedPlugins)
    include_directories("${INTERNAL_SUPERLU_PATH}/SRC/")
	add_library(SuperLU SHARED ${srcSuperLU})	
	target_link_libraries(SuperLU ug4)
	target_link_libraries(SuperLU "${INTERNAL_SUPERLU_PATH}/lib/libsuperlu_4.3.a")
endif(buildEmbeddedPlugins)
