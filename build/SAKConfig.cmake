#####
##
## Copyright (C) 2014 by Luca Heltai
##
#####


##########################################################################
##                                                                      ##
##                The deal.II sak library onfiguration file             ##
##                                                                      ##
##########################################################################


#
# General information
#

SET(SAK_PACKAGE_NAME "")
SET(SAK_PACKAGE_VERSION "")
SET(SAK_PACKAGE_VENDOR "")
SET(SAK_PACKAGE_DESCRIPTION "")

SET(SAK_VERSION_MAJOR "1")
SET(SAK_VERSION_MINOR "0")
SET(SAK_VERSION "1.0")

SET(SAK_PROJECT_CONFIG_NAME "")

SET(SAK_BUILD_TYPE "DebugRelease")
SET(SAK_BUILD_TYPES "")


#
# Information about the project location
#

SET(SAK_CMAKE_MACROS_RELDIR "")
SET(SAK_COMMON_RELDIR "")
SET(SAK_DOCREADME_RELDIR "")
SET(SAK_DOCHTML_RELDIR "")
SET(SAK_EXAMPLES_RELDIR "")
SET(SAK_EXECUTABLE_RELDIR "")
SET(SAK_INCLUDE_RELDIR "")
SET(SAK_LIBRARY_RELDIR "")
SET(SAK_PROJECT_CONFIG_RELDIR "")

#
# Determine SAK_PATH from CMAKE_CURRENT_LIST_DIR:
#

SET(SAK_PATH "${CMAKE_CURRENT_LIST_DIR}")
SET(_path "${SAK_PROJECT_CONFIG_RELDIR}")
WHILE(NOT "${_path}" STREQUAL "")
  GET_FILENAME_COMPONENT(SAK_PATH "${SAK_PATH}" PATH)
  GET_FILENAME_COMPONENT(_path "${_path}" PATH)
ENDWHILE()

#
# Print a message after inclusion of this file:
#

SET(SAK_PROJECT_CONFIG_INCLUDED TRUE)

IF(NOT ${SAK_PACKAGE_NAME}_FIND_QUIETLY)
  MESSAGE(STATUS
    "Using the ${SAK_PACKAGE_NAME}-${SAK_PACKAGE_VERSION} installation found at ${SAK_PATH}"
    )
ENDIF()


#
# Include all convenience macros:
#

FILE(GLOB _macro_files
  "${SAK_PATH}/cmake/macro_*.cmake"
  )
FOREACH(file ${_macro_files})
  IF(NOT ${SAK_PACKAGE_NAME}_FIND_QUIETLY)
    MESSAGE(STATUS "Include macro ${file}")
  ENDIF()
  INCLUDE(${file})
ENDFOREACH()


#
# Compiler and linker configuration
#

SET(SAK_CXX_COMPILER "/Library/deal.II-bundle/openmpi-1.6.5/bin/mpic++")

# used for all targets:
SET(SAK_CXX_FLAGS "-Qunused-arguments -pedantic -fpic -Wall -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Wswitch -Wno-long-long   -Wno-dangling-else -Wno-delete-non-virtual-dtor -Wno-long-long -Wno-newline-eof -Wno-unused-function -Wno-unused-private-field -Wno-unused-variable -Wno-unsupported-friend -std=c++11 -Wno-parentheses -Wno-long-long -Wno-unused -Wno-extra -Wno-overloaded-virtual -Wno-long-long")

# _additionally_ used for debug targets:
SET(SAK_CXX_FLAGS_DEBUG "-O0 -ggdb")

# _additionally_ used for release targets:
SET(SAK_CXX_FLAGS_RELEASE "-O2 -funroll-loops -funroll-all-loops -fstrict-aliasing -Wno-unused")

# used for all targets:
SET(SAK_LINKER_FLAGS " ")

# _additionally_ used for debug targets:
SET(SAK_LINKER_FLAGS_DEBUG "")

# _additionally_ used for release targets:
SET(SAK_LINKER_FLAGS_RELEASE "")

#
# Information about include directories and libraries
#

SET(SAK_INCLUDE_DIRS "")

# Full list of libraries for the debug target:
SET(SAK_LIBRARIES_DEBUG "")

# Full list of libraries for the release target:
SET(SAK_LIBRARIES_RELEASE "")

# Full list of libraries with "debug" and "optimized" keywords for easy use with TARGET_LINK_LIBRARIES:
SET(SAK_LIBRARIES "")


