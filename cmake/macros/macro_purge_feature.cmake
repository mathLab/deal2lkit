## ---------------------------------------------------------------------
##
## Copyright (C) 2014 by the deal2lkit authors
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
# Remove all cached and non cached variables associated with a feature.
#
# Usage:
#     PURGE_FEATURE(feature)
#

MACRO(PURGE_FEATURE _feature)
  #
  # uncached:
  #
  FOREACH(_var ${D2K_LIST_SUFFIXES} ${D2K_STRING_SUFFIXES})
    IF(NOT _var MATCHES BUNDLED)
      SET(${_feature}_${_var})
    ENDIF()
  ENDFOREACH()

  UNSET(${_feature}_FOUND)
  UNSET(${_feature}_VERSION)

  #
  # cached:
  #
  FOREACH(_var ${${_feature}_CLEAR_VARIABLES})
    SET(${_var})
    UNSET(${_var} CACHE)
  ENDFOREACH()

  UNSET(${_feature}_CLEAR_VARIABLES CACHE)

  MARK_AS_ADVANCED(CLEAR ${_feature}_DIR ${_feature}_ARCH)
ENDMACRO()
