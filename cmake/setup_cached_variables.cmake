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
# Set up cached variables (prior to the PROJECT(deal2lkit) call)
#
# This file sets up the following cached Options:
#
# General configuration options:
#
#     D2K_ENABLE_TESTING
#     D2K_COMPONENT_EXAMPLES
#     D2K_COMPILE_EXAMPLES
#     D2K_COMPONENT_DOCUMENTATION
#
# Options regarding compilation and linking:
#
#     CMAKE_BUILD_TYPE
#     BUILD_SHARED_LIBS
#     D2K_CXX_FLAGS                    *)
#     D2K_LINKER_FLAGS                 *)
#
#
#
# *)  May also be set via environment variable (CXXFLAGS, LDFLAGS)
#     (a nonempty cached variable has precedence and will not be
#     overwritten by environment)
#


########################################################################
#                                                                      #
#                    General configuration options:                    #
#                                                                      #
########################################################################


If(D2K_HAVE_DOC_DIRECTORY)
  OPTION(D2K_COMPONENT_DOCUMENTATION
    "Enable configuration, build and installation of the documentation. This adds a COMPONENT \"documentation\" to the build system."
    OFF
    )
ENDIF()

OPTION(D2K_COMPONENT_EXAMPLES
  "Enable configuration and installation of the example steps. This adds a COMPONENT \"examples\" to the build system."
  ON
  )


########################################################################
#                                                                      #
#                       Compilation and linking:                       #
#                                                                      #
########################################################################

#
# Setup CMAKE_BUILD_TYPE:
#

SET(CMAKE_BUILD_TYPE
  "DebugRelease"
  CACHE STRING
  "Choose the type of build, options are: Debug, Release and DebugRelease."
  )

# This is cruel, I know. But it is better to only have a known number of
# options for CMAKE_BUILD_TYPE...
IF( NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release" AND
    NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" AND
    NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease" )
  MESSAGE(FATAL_ERROR
    "CMAKE_BUILD_TYPE does neither match Release, Debug, nor DebugRelease!"
    )
ENDIF()

#
# Configuration behaviour:
#


OPTION(D2K_SETUP_DEFAULT_COMPILER_FLAGS
  "Configure sensible default CFLAGS and CXXFLAGS depending on platform, compiler and build target."
  ON
  )
MARK_AS_ADVANCED(D2K_SETUP_DEFAULT_COMPILER_FLAGS)


SET(BUILD_SHARED_LIBS "ON" CACHE BOOL
  "Build a shared library"
  )


SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH "ON" CACHE BOOL
  "Set the rpath of the library to the external link pathes on installation"
  )
MARK_AS_ADVANCED(CMAKE_INSTALL_RPATH_USE_LINK_PATH)


#
# Translate CMake specific variables to deal.II naming:
#

FOREACH(_flag CXX_FLAGS CXX_FLAGS_RELEASE CXX_FLAGS_DEBUG)
  IF(NOT "${CMAKE_${_flag}}" STREQUAL "")
    MESSAGE(STATUS
      "Prepending \${CMAKE_${_flag}} to \${D2K_${_flag}}"
      )
    SET(D2K_${_flag} "${CMAKE_${_flag}} ${D2K_${_flag}}")
  ENDIF()
ENDFOREACH()

FOREACH(_flag LINKER_FLAGS LINKER_FLAGS_DEBUG LINKER_FLAGS_RELEASE)
  IF(NOT "${CMAKE_SHARED_${_flag}}" STREQUAL "")
    MESSAGE(STATUS
      "Prepending \${CMAKE_SHARED_${_flag}} to \${D2K_${_flag}}"
      )
    SET(D2K_${_flag} "${CMAKE_${_flag}} ${D2K_${_flag}}")
  ENDIF()
ENDFOREACH()

#
# Hide all unused CMake variables:
#

SET(D2K_REMOVED_FLAGS
  CMAKE_CXX_FLAGS
  CMAKE_CXX_FLAGS_RELEASE
  CMAKE_CXX_FLAGS_DEBUG
  CMAKE_CXX_FLAGS_MINSIZEREL
  CMAKE_CXX_FLAGS_RELWITHDEBINFO
  CMAKE_C_FLAGS
  CMAKE_C_FLAGS_RELEASE
  CMAKE_C_FLAGS_DEBUG
  CMAKE_C_FLAGS_MINSIZEREL
  CMAKE_C_FLAGS_RELWITHDEBINFO
  CMAKE_Fortran_FLAGS
  CMAKE_Fortran_FLAGS_RELEASE
  CMAKE_Fortran_FLAGS_DEBUG
  CMAKE_Fortran_FLAGS_MINSIZEREL
  CMAKE_Fortran_FLAGS_RELWITHDEBINFO
  CMAKE_SHARED_LINKER_FLAGS
  CMAKE_SHARED_LINKER_FLAGS_DEBUG
  CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL
  CMAKE_SHARED_LINKER_FLAGS_RELEASE
  CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO
  )
FOREACH(_flag ${D2K_REMOVED_FLAGS})
  # Go away...
  SET(${_flag} ${${_flag}} CACHE INTERNAL "" FORCE)
  # Also set it to an empty string for the configuration run so that it
  # does not confuse the build system (to unset is not an option - it is
  # cached...)
  SET(${_flag} "")
ENDFOREACH()

#
# Promote our configuration variables to cache:
#

SET(D2K_USED_FLAGS
  D2K_CXX_FLAGS
  D2K_CXX_FLAGS_DEBUG
  D2K_CXX_FLAGS_RELEASE
  D2K_LINKER_FLAGS
  D2K_LINKER_FLAGS_DEBUG
  D2K_LINKER_FLAGS_RELEASE
  )
FOREACH(_flag ${D2K_USED_FLAGS})
  #
  # Promote to cache:
  #
  SET(${_flag} "${${_flag}}" CACHE STRING
    "The user supplied cache variable will be appended _at the end_ of the configuration step to the auto generated ${_flag} variable"
    )
  MARK_AS_ADVANCED(${_flag})

  #
  # The order of compiler and linker flags is important. In order to
  # provide an override mechanism we have to save the initial (cached)
  # variable at this point and clear it.
  # ${flags}_SAVED will be appended to ${flags} again in
  # setup_finalize.cmake (called at the end of the main CMakeLists.txt
  # file).
  #
  SET(${_flag}_SAVED ${${_flag}})
  SET(${_flag} "")
ENDFOREACH()

FOREACH(_variable
  D2K_DEFINITIONS
  D2K_DEFINITIONS_DEBUG
  D2K_DEFINITIONS_RELEASE
  )
  #
  # Promote to cache:
  #
  SET(${_variable} ${${_variable}} CACHE STRING
    "Additional, user supplied compile definitions"
    )
  MARK_AS_ADVANCED(${_variable})
ENDFOREACH()


#
# Finally, read in CXXFLAGS and LDFLAGS from environment and prepend them
# to the saved variables:
#
# Also strip leading and trailing whitespace from linker flags to make
# old cmake versions happy
#
SET(D2K_CXX_FLAGS_SAVED "$ENV{CXXFLAGS} ${D2K_CXX_FLAGS_SAVED}")
STRING(STRIP "${D2K_CXX_FLAGS_SAVED}" D2K_CXX_FLAGS_SAVED)
SET(D2K_LINKER_FLAGS_SAVED "$ENV{LDFLAGS} ${D2K_LINKER_FLAGS_SAVED}")
STRING(STRIP "${D2K_LINKER_FLAGS_SAVED}" D2K_LINKER_FLAGS_SAVED)
UNSET(ENV{CXXFLAGS})
UNSET(ENV{LDFLAGS})


########################################################################
#                                                                      #
#                Components and miscellaneous setup:                   #
#                                                                      #
########################################################################


OPTION(D2K_DOXYGEN_USE_MATHJAX
  "If set to ON, doxygen documentation is generated using mathjax"
  OFF
  )
MARK_AS_ADVANCED(D2K_DOXYGEN_USE_MATHJAX)



########################################################################
#                                                                      #
#                               Finalize:                              #
#                                                                      #
########################################################################

#
# We do not support installation into the binary directory any more ("too
# much pain, not enough profit"):
#

IF("${CMAKE_BINARY_DIR}" STREQUAL "${CMAKE_INSTALL_PREFIX}")
  MESSAGE(FATAL_ERROR "
Error CMAKE_INSTALL_PREFIX is equal to CMAKE_BINARY_DIR.
It is not possible to install into the build directory. Please set
CMAKE_INSTALL_PREFIX to a designated install directory different than
CMAKE_BINARY_DIR.
(Please note that you can use deal.II directly out of a build directory
without the need to install it, if this is what you tried to do.)
"
    )
ENDIF()

#
# Compatibility renaming:
#

IF(DEFINED D2K_HAVE_CXX11_FLAG AND NOT D2K_HAVE_CXX11_FLAG)
  SET(D2K_WITH_CXX11 FALSE CACHE BOOL "" FORCE)
ENDIF()

#
# Miscellaneous renaming:
#

GET_CMAKE_PROPERTY(_res VARIABLES)
FOREACH(_var ${_res})
  #
  # Rename (ALLOW|WITH|FORCE|COMPONENT)_* by D2K_(ALLOW|WITH|FORCE|COMPONENT)_*
  #
  FOREACH(_match ALLOW_ WITH_ FORCE_ COMPONENT_)
    IF(_var MATCHES "^${_match}")
      SET(D2K_${_var} ${${_var}} CACHE BOOL "" FORCE)
      UNSET(${_var} CACHE)
    ENDIF()
  ENDFOREACH()

  #
  # Same for components:
  #
  IF(_var MATCHES "^(DOCUMENTATION|EXAMPLES|PACKAGE|PARAMETER_GUI)")
    SET(D2K_COMPONENT_${_var} ${${_var}} CACHE BOOL "" FORCE)
    UNSET(${_var} CACHE)
  ENDIF()

ENDFOREACH()
