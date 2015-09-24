## ---------------------------------------------------------------------
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
# A small wrapper around FIND_LIBRARY to be a bit more verbose
#

MACRO(D2K_FIND_LIBRARY _library_name)
  FIND_LIBRARY(${_library_name} ${ARGN})

  IF(${_library_name} MATCHES "-NOTFOUND")
    MESSAGE(STATUS "${_library_name} not found! Call:")
    TO_STRING(_str ${ARGN})
    MESSAGE(STATUS "    FIND_LIBRARY(${_library_name} ${_str})")
  ELSE()
    MESSAGE(STATUS "Found ${_library_name}")
  ENDIF()
ENDMACRO()
