#####
##
## Copyright (C) 2013 Luca Heltai
##
## This file is part of the deal.ii SAK library.
##
## Author: Luca Heltai <luca.heltai@sissa.it>
##	   from original design by 
##	   Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
MESSAGE(STATUS "Building deal.II SAK library.")
PROJECT(SAK)

set(SAK_VERSION_MAJOR 1)
set(SAK_VERSION_MINOR 0)
set(SAK_VERSION
  ${SAK_VERSION_MAJOR}.${MATHLAB_VERSION_MINOR})

SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

# Include all macros. Some of them are copied brutally from deal.II
FILE(GLOB _cmake_files "cmake/macros/*.cmake")
FOREACH(_file ${_cmake_files})
  INCLUDE(${_file})
ENDFOREACH()

# These ones are for GSL, OCE and SUNDIALS
LIST(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
  ${CMAKE_SOURCE_DIR}/cmake/modules/
  )

find_package(deal.II REQUIRED HINTS $ENV{DEAL_II_DIR}) 
# find_package(OpenCASCADE MODULE REQUIRED) 
# find_package(GSL MODULE REQUIRED) 
# find_package(SUNDIALS MODULE REQUIRED) 

include_directories(
    ${DEAL_II_INCLUDE_DIRS}
#    ${OpenCASCADE_INCLUDE_DIRS}
#    ${SUNDIALS_INCLUDE_DIRS}
#    ${GSL_INCLUDE_DIRS}
    include 
) 

set(CMAKE_BUILD_TYPE "DebugRelease")

SET_IF_EMPTY(_install $ENV{SAK_SRC})
SET_IF_EMPTY(_install "/usr/local/sak")

SET(CMAKE_INSTALL_PREFIX ${_install})

SET(CMAKE_CXX_FLAGS ${DEAL_II_CXX_FLAGS})
SET(CMAKE_CXX_FLAGS_RELEASE ${DEAL_II_CXX_FLAGS_RELEASE})
SET(CMAKE_CXX_FLAGS_DEBUG ${DEAL_II_CXX_FLAGS_DEBUG})

file(GLOB cc-files source/*.cc) 
file(GLOB h-files  include/*h)

# Almost verbatim identical to DEAL_II_ADD_LIBRARY. I use this to 
# set the properties of both dealii-sak.debug and dealii-sak.release targets, 
# 

SAK_ADD_LIBRARY(dealii-sak SHARED ${cc-files})

install(DIRECTORY ${_install}/include/ DESTINATION include
        FILES_MATCHING PATTERN "*.h")


#############################################################
# Now try a SAKConfig.cmake
#############################################################
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# Create the SAKConfig.cmake and MATHLABConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR 
    "${INSTALL_CMAKE_DIR}"
    "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS 
    "${PROJECT_SOURCE_DIR}" 
    "${PROJECT_BINARY_DIR}")

configure_file(cmake/SAKConfig.cmake.in
  "${PROJECT_BINARY_DIR}/SAKConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${SAK_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(cmake/SAKConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/SAKConfig.cmake" @ONLY)
# ... for both
configure_file(cmake/SAKConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/SAKConfigVersion.cmake" @ONLY)
 
# Install the SAKConfig.cmake and MATHLABConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/SAKConfig.cmake"
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/SAKConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
