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
# This file implements the D2K_SETUP_TARGET macro, which is
# part of the deal.II library.
#
# Usage:
#       D2K_SETUP_TARGET(target)
#       D2K_SETUP_TARGET(target DEBUG|RELEASE)
#
# This appends necessary include directories, linker flags, compile flags
# and compile definitions and the deal.II library link interface to the
# given target. In particular:
#
# INCLUDE_DIRECTORIES is appended with
#   "${D2K_INCLUDE_DIRS}"
#
# COMPILE_FLAGS is appended with
#   "${D2K_CXX_FLAGS} ${D2K_CXX_FLAGS_<build type>}"
#
# LINK_FLAGS is appended with
#   "${D2K_LINKER_FLAGS ${D2K_LINKER_FLAGS_<build type>}"
#
# COMPILE_DEFINITIONS is appended with
#   "${D2K_USER_DEFINITIONS};${D2K_USER_DEFINITIONS_<build type>}"
#
# If no "DEBUG" or "RELEASE" keyword is specified after the target, the
# current CMAKE_BUILD_TYPE determines which compiler and linker flags as
# well as compile definitions to use and against which deal.II library it
# should be linked against.
#
# If the requested build type is not available (e.g. DEBUG request but
# deal.II was compiled with release mode only), the other available will be
# used instead.
#

MACRO(D2K_SETUP_TARGET _target)

  IF(NOT D2K_PROJECT_CONFIG_INCLUDED)
    MESSAGE(FATAL_ERROR
      "\nD2K_SETUP_TARGET can only be called in external projects after "
      "the inclusion of deal2lkitConfig.cmake. It is not intended for "
      "internal use.\n\n"
      )
  ENDIF()

  IF(NOT D2K_TARGET_CONFIG_INCLUDED)
    INCLUDE(${D2K_TARGET_CONFIG})
    SET(D2K_TARGET_CONFIG_INCLUDED TRUE)
  ENDIF()

  # Necessary for setting INCLUDE_DIRECTORIES via SET_PROPERTY
  CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

  #
  # Every build type that (case insensitively) matches "debug" is
  # considered a debug build:
  #
  SET(_build "RELEASE")
  STRING(TOLOWER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
  IF("${_cmake_build_type}" MATCHES "debug")
    SET(_build "DEBUG")
  ENDIF()

  #
  # Override _on_debug_build if ${ARGN} is set:
  #
  IF("${ARGN}" MATCHES "^(DEBUG|RELEASE)$")
    SET(_build "${ARGN}")
  ENDIF()

  #
  # We can only append DEBUG link flags and compile definitions if deal.II
  # was built with the Debug or DebugRelease build type. So test for this:
  #
  IF("${_build}" STREQUAL "DEBUG" AND NOT D2K_BUILD_TYPE MATCHES "Debug")
    SET(_build "RELEASE")
  ENDIF()

  SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
    INCLUDE_DIRECTORIES "${D2K_INCLUDE_DIRS}"
    )

  #
  # Set up the link interface:
  #
  GET_PROPERTY(_type TARGET ${_target} PROPERTY TYPE)
  IF(NOT "${_type}" STREQUAL "OBJECT_LIBRARY")
    TARGET_LINK_LIBRARIES(${_target} ${D2K_TARGET_${_build}})
  ENDIF()

  #
  # If D2K_STATIC_EXECUTABLE is set, switch the final link type to
  # static:
  #
  IF(D2K_STATIC_EXECUTABLE)
    SET_PROPERTY(TARGET ${_target} PROPERTY
      LINK_SEARCH_END_STATIC TRUE
      )
  ENDIF()

  DEAL_II_SETUP_TARGET(${_target} ${ARGN})

ENDMACRO()
