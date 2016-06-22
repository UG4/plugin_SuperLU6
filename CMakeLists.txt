# Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
# Author: Martin Rupp
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

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

set(pluginName	SuperLU)
set(SOURCES 	super_lu_bridge.cpp super_lu.cpp super_lu.h)



cmake_minimum_required(VERSION 2.6)
project(UG_PLUGIN_${pluginName})
# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)


# detect path for SuperLU-Library
if(SUPERLU_PATH)
    # the user can supply a SUPER_LU path for his own compiled SuperLU   
    set(INTERNAL_SUPERLU_PATH "${SUPERLU_PATH}")
    message(STATUS "Info: SuperLU is ON, using user provided path SUPERLU_PATH=${INTERNAL_SUPERLU_PATH}")
else(SUPERLU_PATH)
    UGInstall("SuperLU" INTERNAL_SUPERLU_PATH)
    message(STATUS "Info: SuperLU is ON, using uginstall provided path ${INTERNAL_SUPERLU_PATH} (Specify your own path by -DSUPERLU_PATH=... if necessary.)")    
endif(SUPERLU_PATH)


# clear library variable to make sure it is updated if necessary
SET(SLU_LIB SLU_LIB-NOTFOUND)

# we only support version 4.3
#find_library(SLU_LIB NAMES superlu_4.3 superlu_4.2 superlu_4.1 superlu_4.0 superlu_5.1.1 superlu_5.1 

find_library(SLU_LIB NAMES superlu_4.3
			 PATHS ${INTERNAL_SUPERLU_PATH} PATH_SUFFIXES lib
			 NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

if ("${SLU_LIB}" STREQUAL "SLU_LIB-NOTFOUND")
	message(FATAL_ERROR "SuperLU library was not found in ${INTERNAL_SUPERLU_PATH}")
endif()


# Add to build process 
if(buildEmbeddedPlugins)
	ExportSources(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
    ExportDependencies("${SLU_LIB}")
    if(EXISTS "${INTERNAL_SUPERLU_PATH}/lib/blas.a")
        ExportDependencies("${INTERNAL_SUPERLU_PATH}/lib/blas.a")
    endif(EXISTS "${INTERNAL_SUPERLU_PATH}/lib/blas.a")
    ExportIncludes("${INTERNAL_SUPERLU_PATH}/SRC/")
else(buildEmbeddedPlugins)
    include_directories("${INTERNAL_SUPERLU_PATH}/SRC/")
	add_library(SuperLU SHARED ${SOURCES})	
	target_link_libraries(SuperLU ug4)
	target_link_libraries(SuperLU "${SLU_LIB}")
endif(buildEmbeddedPlugins)
