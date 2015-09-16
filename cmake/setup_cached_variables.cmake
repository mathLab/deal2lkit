## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2015 by the deal2lkit authors
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

OPTION(D2K_ENABLE_TESTING 
  "Enable configuration and installation of the tests." 
	OFF
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



SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH "ON" CACHE BOOL
  "Set the rpath of the library to the external link pathes on installation"
  )
MARK_AS_ADVANCED(CMAKE_INSTALL_RPATH_USE_LINK_PATH)


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
(Please note that you can use deal2lkit directly out of a build directory
without the need to install it, if this is what you tried to do.)
"
    )
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
