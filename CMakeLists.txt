## ---------------------------------------------------------------------
##
## This file is copyrighted by the deal2lkit authors and by the deal.II 
## authors (see below for the original copyright in the deal.II library.)
## 
## The structure of the cmake files are the same of those of the 
## deal.II library and they have been modified in order to make them
## compatible with the deal2lkit library.
##
## This file is part of the deal2lkit library.
##
## The deal2lkit library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal2lkit distribution.
## ---------------------------------------------------------------------

## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

##                                                                    ##
#            The cmake build system for the deal2lkit project          #
##                                                                    ##

########################################################################
#                                                                      #
#                            Configuration:                            #
#                                                                      #
########################################################################

#
# General configuration for cmake:
#
MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

#
# We support all policy changes up to version 3.1.0. Thus, explicitly set
# all policies CMP0001 - CMP0054 to new for version 3.1 (and later) to
# avoid some unnecessary warnings.
#
IF( "${CMAKE_VERSION}" VERSION_EQUAL "3.1" OR
    "${CMAKE_VERSION}" VERSION_GREATER "3.1" )
  CMAKE_POLICY(VERSION 3.1.0)
ENDIF()

IF(POLICY CMP0037)
  # allow to override "test" target for quick tests
  CMAKE_POLICY(SET CMP0037 OLD)
ENDIF()

FIND_PACKAGE(deal.II 8.3 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )

DEAL_II_INITIALIZE_CACHED_VARIABLES()
DEAL_II_QUERY_GIT_INFORMATION()

# deal2lkit requires c++11 support
IF(NOT DEAL_II_WITH_CXX11)
  MESSAGE(FATAL_ERROR "
  Error! The deal.II library found at ${DEAL_II_PATH} was not configured with
  DEAL_II_WITH_CXX11 = ON
  Please recompile deal.II enabling this feature."
  )
ENDIF()

#
# The INCLUDE() and FIND_PACKAGE() search first in these directories
#
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)


#
# Load all macros:
#
FILE(GLOB _macro_files "cmake/macros/*.cmake")
MESSAGE(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/load_cmake_macros.cmake")
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/load_cmake_macros.cmake)
FOREACH(_file ${_macro_files})
  MESSAGE(STATUS "Include ${_file}")
  INCLUDE(${_file})
ENDFOREACH()

#D2K_QUERY_GIT_INFORMATION()

#
# Check for the existence of various optional folders:
#
IF(EXISTS ${CMAKE_SOURCE_DIR}/doc/CMakeLists.txt)
  SET(D2K_HAVE_DOC_DIRECTORY TRUE)
ENDIF()

IF(EXISTS ${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt)
  SET(D2K_HAVE_TESTS_DIRECTORY TRUE)
ENDIF()

#
# We have to initialize some cached variables before PROJECT is called, so
# do it at this point:
#
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_cached_variables.cmake)

#
# Now, set the project and set up the rest:
#
PROJECT(deal2lkit CXX C)

INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_deal2lkit.cmake)


#
# Feature configuration:
#
FILE(GLOB _configure_files "cmake/configure/configure_*.cmake")
LIST(SORT _configure_files) # make sure to include in alphabetical order
FOREACH(_file ${_configure_files})
  INCLUDE(${_file})
ENDFOREACH()


#
# Finalize the configuration:
#
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_custom_targets.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_finalize.cmake)
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/setup_write_config.cmake)

########################################################################
#                                                                      #
#                     Compilation and installation:                    #
#                                                                      #
########################################################################

MESSAGE(STATUS "")
MESSAGE(STATUS "Configuring done. Proceed to target definitions now.")

ADD_SUBDIRECTORY(include)


#ADD_SUBDIRECTORY(source) 

INCLUDE_DIRECTORIES(
  ${CMAKE_BINARY_DIR}/include
  ${CMAKE_SOURCE_DIR}/include
  ${D2K_INCLUDE_DIRS}
  )
file(GLOB _library_source_files source/*.cc include/deal2lkit/*h)

FOREACH(_BUILD_TYPE ${D2K_BUILD_TYPES})
  SET(_lib ${D2K_BASE_NAME}${D2K_${_BUILD_TYPE}_SUFFIX})
  ADD_LIBRARY(${_lib}
    SHARED
    ${_library_source_files}
    )


  # link external libraries
  TARGET_LINK_LIBRARIES(${_lib}
    ${D2K_LIBRARIES_${_BUILD_TYPE}}
    ${D2K_LIBRARIES}
    )
  SET(TEST_LIBRARIES_${_BUILD_TYPE} ${_lib})

  SET_TARGET_PROPERTIES(${_lib} PROPERTIES MACOSX_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${D2K_LIBRARY_RELDIR}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${D2K_LIBRARY_RELDIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${D2K_LIBRARY_RELDIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${D2K_EXECUTABLE_RELDIR}"
   )

  DEAL_II_SETUP_TARGET(${_lib} ${_BUILD_TYPE})


  FILE(MAKE_DIRECTORY
    ${CMAKE_BINARY_DIR}/${D2K_PROJECT_CONFIG_RELDIR}
    )
  EXPORT(TARGETS ${_lib}
    FILE ${CMAKE_BINARY_DIR}/${D2K_PROJECT_CONFIG_RELDIR}/${D2K_PROJECT_CONFIG_NAME}Targets.cmake
    APPEND
    )


  INSTALL(TARGETS ${D2K_BASE_NAME}${D2K_${_BUILD_TYPE}_SUFFIX}
    COMPONENT library
    EXPORT ${D2K_PROJECT_CONFIG_NAME}Targets
    RUNTIME DESTINATION ${D2K_EXECUTABLE_RELDIR}
    LIBRARY DESTINATION ${D2K_LIBRARY_RELDIR}
    ARCHIVE DESTINATION ${D2K_LIBRARY_RELDIR}
    )
ENDFOREACH()
INSTALL(EXPORT ${D2K_PROJECT_CONFIG_NAME}Targets
  DESTINATION ${D2K_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )

MESSAGE(STATUS "Setting up library - Done")



ADD_SUBDIRECTORY(cmake/config) # has to be included after source
ADD_SUBDIRECTORY(examples)

#Do not automatically include the testsuite subproject on Windows.
IF(D2K_HAVE_TESTS_DIRECTORY AND NOT DEAL_II_MSVC AND D2K_ENABLE_TESTING)
  ADD_SUBDIRECTORY(tests)
#  ENABLE_TESTING()
ENDIF()

IF(D2K_COMPONENT_DOCUMENTATION)
  ADD_SUBDIRECTORY(doxygen)
ENDIF()

#
# And finally, print the configuration:
#
FILE(READ ${CMAKE_BINARY_DIR}/summary.log D2K_LOG_SUMMARY)
MESSAGE("${D2K_LOG_SUMMARY}")
