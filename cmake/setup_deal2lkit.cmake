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
# Set up deal2lkit specific definitions
#
# Definitions marked with *) can be overriden by defining them to cache
# prior to the call of this file. This is done with the help of the
# SET_IF_EMPTY macro.
#
# General information about deal2lkit:
#
#     D2K_PACKAGE_NAME            *)
#     D2K_PACKAGE_VERSION         *)
#     D2K_PACKAGE_VENDOR          *)
#     D2K_PACKAGE_DESCRIPTION     *)
#     D2K_VERSION_MAJOR
#     D2K_VERSION_MINOR
#     D2K_VERSION_SUBMINOR
#     D2K_VERSION
#
# Information about paths, install locations and names:
#
#     D2K_PROJECT_CONFIG_NAME     *)
#     D2K_BASE_NAME               *)
#     D2K_DEBUG_SUFFIX            *)
#     D2K_RELEASE_SUFFIX          *)
#
#     D2K_EXECUTABLE_RELDIR       *)
#     D2K_INCLUDE_RELDIR          *)
#     D2K_LIBRARY_RELDIR          *)
#     D2K_PROJECT_CONFIG_RELDIR   *)
#     D2K_SHARE_RELDIR            *)
#     D2K_DOCREADME_RELDIR        *)
#     D2K_DOCHTML_RELDIR          *)
#     D2K_EXAMPLES_RELDIR         *)
#
#     D2K_BUILD_TYPES
#
# *)  Can be overwritten by the command line via -D<...>
#

########################################################################
#                                                                      #
#                  General information about deal2lkit:                #
#                                                                      #
########################################################################

SET_IF_EMPTY(D2K_PACKAGE_NAME "deal2lkit")

SET_IF_EMPTY(D2K_PACKAGE_VENDOR
  "The deal2lkit Authors <https://github.com/mathLab/deal2lkit/>"
  )
SET_IF_EMPTY(D2K_PACKAGE_DESCRIPTION
	"A Toolkit library for deal.II <http://www.dealii.org/>"
  )

FILE(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" _version LIMIT_COUNT 1)
SET_IF_EMPTY(D2K_PACKAGE_VERSION "${_version}")

#
# We expect a version number of the form "X.Y.Z", where X and Y are always
# numbers and Z is either a third number (for a release version) or a short
# string.
#
STRING(REGEX REPLACE "^([0-9]+)\\..*" "\\1"
  D2K_VERSION_MAJOR "${D2K_PACKAGE_VERSION}"
  )
STRING(REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1"
  D2K_VERSION_MINOR "${D2K_PACKAGE_VERSION}"
  )

#
# If Z is not a number, replace it with "0", otherwise extract version
# number:
#
IF(D2K_PACKAGE_VERSION MATCHES "^[0-9]+\\.[0-9]+.*\\.[0-9]+.*")
  STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+.*\\.([0-9]+).*" "\\1"
    D2K_VERSION_SUBMINOR "${D2K_PACKAGE_VERSION}"
    )
ELSE()
  SET(D2K_VERSION_SUBMINOR "0")
ENDIF()
SET(D2K_VERSION ${D2K_VERSION_MAJOR}.${D2K_VERSION_MINOR}.${D2K_VERSION_SUBMINOR})


########################################################################
#                                                                      #
#         Information about paths, install locations and names:        #
#                                                                      #
########################################################################

SET(D2K_PROJECT_CONFIG_NAME "${D2K_PACKAGE_NAME}")

STRING(REPLACE "." "_" _base_name "${D2K_PACKAGE_NAME}")
SET_IF_EMPTY(D2K_BASE_NAME "${_base_name}")
SET_IF_EMPTY(D2K_DEBUG_SUFFIX ".g")
SET_IF_EMPTY(D2K_RELEASE_SUFFIX "")

#
# Try to obey the FSHS as close as possible ...
#
SET_IF_EMPTY(D2K_EXECUTABLE_RELDIR "bin")
SET_IF_EMPTY(D2K_INCLUDE_RELDIR "include")
SET_IF_EMPTY(D2K_LIBRARY_RELDIR "lib${LIB_SUFFIX}")
SET_IF_EMPTY(D2K_PROJECT_CONFIG_RELDIR "${D2K_LIBRARY_RELDIR}/cmake/${D2K_PROJECT_CONFIG_NAME}")
SET_IF_EMPTY(D2K_SHARE_RELDIR "share/${D2K_PACKAGE_NAME}")
#
# ... but install the documentation into prominent places:
#
SET_IF_EMPTY(D2K_DOCREADME_RELDIR "./")
SET_IF_EMPTY(D2K_DOCHTML_RELDIR "doc")
SET_IF_EMPTY(D2K_EXAMPLES_RELDIR "examples")

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  LIST(APPEND D2K_BUILD_TYPES "DEBUG")
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES "Release")
  LIST(APPEND D2K_BUILD_TYPES "RELEASE")
ENDIF()


SET(D2K_LIST_SUFFIXES
  DEFINITIONS DEFINITIONS_RELEASE DEFINITIONS_DEBUG
  USER_DEFINITIONS USER_DEFINITIONS_RELEASE USER_DEFINITIONS_DEBUG
  INCLUDE_DIRS USER_INCLUDE_DIRS BUNDLED_INCLUDE_DIRS
  LIBRARIES LIBRARIES_RELEASE LIBRARIES_DEBUG
  )

SET(D2K_STRING_SUFFIXES
  CXX_FLAGS CXX_FLAGS_RELEASE CXX_FLAGS_DEBUG
  LINKER_FLAGS LINKER_FLAGS_RELEASE LINKER_FLAGS_DEBUG
  )


