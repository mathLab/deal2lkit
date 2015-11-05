##-----------------------------------------------------------
##
##    Copyright (C) 2015 by the deal2lkit authors
##
##    This file is part of the deal2lkit library.
##
##    The deal2lkit library is free software; you can use it, redistribute
##    it, and/or modify it under the terms of the GNU Lesser General
##    Public License as published by the Free Software Foundation; either
##    version 2.1 of the License, or (at your option) any later version.
##    The full text of the license can be found in the file LICENSE at
##    the top level of the deal2lkit distribution.
##
##-----------------------------------------------------------

##
# CMake script for small project
#
# If you set the environemnt variable DEAL2LKIT_DIR or D2K_DIR, 
# everything will work out of the box
##

# Set the name of the project and target
# If your application follows the structure above, you don't need to 
# specify anything else. 
SET(TARGET poisson)

############################################################
# Normally you shouldn't need to change anything below.
############################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.4 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES()


PROJECT(${TARGET})
FIND_PACKAGE(deal2lkit 1.0 REQUIRED
  HINTS ${D2K_DIR} $ENV{D2K_DIR} $ENV{DEAL2LKIT_DIR}
  )
D2K_INITIALIZE_CACHED_VARIABLES()

# We add one library and one target for each type of deal.II library
# we found. If you compiled deal.II with both Release and Debug
# mode, this will generate both Release and Debug programs for you.
# The debug library and program are postfixed with ".g"
SET(_d2_build_types "Release;Debug")
SET(Release_postfix "")
SET(Debug_postfix ".g")

FOREACH(_build_type ${_d2_build_types})
  # Postfix to use everywhere
  SET(_p "${${_build_type}_postfix}")
  # Only build this type, if deal.II was compiled with it.
  IF(CMAKE_BUILD_TYPE MATCHES "${_build_type}" AND
      DEAL_II_BUILD_TYPE MATCHES "${_build_type}"  AND
      D2K_BUILD_TYPE MATCHES "${_build_type}")

    MESSAGE("-- Found ${_build_type} version of deal.II.")
    MESSAGE("-- Found ${_build_type} version of deal2lkit.")

    SET(_exe "${TARGET}${${_build_type}_postfix}")
    MESSAGE("-- Configuring executable ${_exe}")
    ADD_EXECUTABLE(${_exe} ${TARGET}.cc)
    D2K_SETUP_TARGET(${_exe} ${_BUILD_TYPE})
  ENDIF()
ENDFOREACH()
