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
##
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

#
# This file implements the D2K_INITIALIZE_VARIABLES macro, which is
# part of the deal2lkit library.
#
# Usage:
#       D2K_INITIALIZE_CACHED_VARIABLES()
#
# This sets some cached variables to the values used for compiling the
# deal2lkit library.
#
# This macro has to be called before PROJECT()!
#

MACRO(D2K_INITIALIZE_CACHED_VARIABLES)

  IF(NOT D2K_PROJECT_CONFIG_INCLUDED)
    MESSAGE(FATAL_ERROR
      "\nD2K_INITIALIZE_CACHED_VARIABLES can only be called in external "
      "projects after the inclusion of deal2lkitConfig.cmake. It is not "
      "intended for internal use.\n\n"
      )
  ENDIF()

  #
  # Set build type according to available libraries
  #
  IF(D2K_BUILD_TYPE MATCHES "Debug")
    SET(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug, Release"
      )
  ELSE()
    SET(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug, Release"
      )
  ENDIF()

  #
  # Reset build type if unsupported, i.e. if it is not (case insensitively
  # equal to Debug or Release or unsupported by the current build type:
  #
  STRING(TOLOWER "${CMAKE_BUILD_TYPE}" _cmake_build_type)

  IF(NOT "${_cmake_build_type}" MATCHES "^(debug|release|debugrelease)$")

    IF("${D2K_BUILD_TYPE}" STREQUAL "DebugRelease")
      SET(_new_build_type "Debug")
    ELSE()
      SET(_new_build_type "${D2K_BUILD_TYPE}")
    ENDIF()

    MESSAGE(
      "###
      #
      #  WARNING:
      #
      #  CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\" unsupported by current installation!
      #  deal2lkit was built with CMAKE_BUILD_TYPE \"${D2K_BUILD_TYPE}\".
      #
      #  CMAKE_BUILD_TYPE is forced to \"${_new_build_type}\".
      #
      ###"
      )
    SET(CMAKE_BUILD_TYPE "${_new_build_type}" CACHE STRING
      "Choose the type of build, options are: Debug, Release"
      FORCE
      )

  ENDIF()


  SET(CMAKE_CXX_COMPILER ${D2K_CXX_COMPILER} CACHE STRING
    "CXX Compiler.")

  SET(CMAKE_CXX_FLAGS "" CACHE STRING
    "Flags used by the compiler during all build types."
    )

  SET(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING
    "Flags used by the compiler during debug builds."
    )

  SET(CMAKE_CXX_FLAGS_RELEASE "" CACHE STRING
    "Flags used by the compiler during release builds."
    )

  MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)

ENDMACRO()

