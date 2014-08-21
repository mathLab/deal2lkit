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

include_directories(
    ${DEAL_II_INCLUDE_DIRS}
    include 
) 

set(CMAKE_BUILD_TYPE "DebugRelease")

SET(CMAKE_CXX_FLAGS ${DEAL_II_CXX_FLAGS})
SET(CMAKE_CXX_FLAGS_RELEASE ${DEAL_II_CXX_FLAGS_RELEASE})
SET(CMAKE_CXX_FLAGS_DEBUG ${DEAL_II_CXX_FLAGS_DEBUG})

file(GLOB cc-files source/*.cc) 
file(GLOB h-files  include/*h)

# Almost verbatim identical to DEAL_II_ADD_LIBRARY. I use this to 
# set the properties of both dealii-sak.debug and dealii-sak.release targets, 
# 

SAK_ADD_LIBRARY(dealii-sak SHARED ${cc-files})

install(DIRECTORY include/ DESTINATION include
        FILES_MATCHING PATTERN "*.h")


#-- Add an Option to toggle the generation of the API documentation
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
if(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)
  if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR 
      "Doxygen is needed to build the documentation. Please install it correctly or set BUILD_DOCUMENTATION to OFF")
  endif()
  #-- Configure the Template Doxyfile for our specific project
  configure_file(Doxyfile.in 
                 ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
  #-- Add a custom target to run Doxygen when ever the project is built
  add_custom_target (docs 
    COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
    SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
endif()
