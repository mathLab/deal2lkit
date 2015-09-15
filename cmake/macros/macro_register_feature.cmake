## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2014 by the deal.II authors
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
# This macro is used for the feature configuration in deal.II. It adds
# individual FEATURE_* configuration variables to the corresponding
# D2K_* variables
#
# Usage:
#     REGISTER_FEATURE(feature)
#
# This macro will add
#
#   <FEATURE>_LIBRARIES (respecting general, optimized, debug keyword)
#
# and all other suffixes defined in D2K_LIST_SUFFIXES and
# D2K_STRING_SUFFIXES to the corresponding D2K_* variables
#

MACRO(REGISTER_FEATURE _feature)

  IF(DEFINED ${_feature}_LIBRARIES)
    #
    # Add ${_feature}_LIBRARIES to
    #   D2K_LIBRARIES
    #   D2K_LIBRARIES_DEBUG
    #   D2K_LIBRARIES_RELEASE
    # depending on the "optmized", "debug" or "general" keyword
    #
    SET(_toggle "general")
    FOREACH(_tmp ${${_feature}_LIBRARIES})
      IF( "${_tmp}" STREQUAL "debug" OR
          "${_tmp}" STREQUAL "optimized" OR
          "${_tmp}" STREQUAL "general" )
        SET(_toggle "${_tmp}")
      ELSE()
        IF("${_toggle}" STREQUAL "general")
          LIST(APPEND D2K_LIBRARIES ${_tmp})
        ELSEIF("${_toggle}" STREQUAL "debug")
          LIST(APPEND D2K_LIBRARIES_DEBUG ${_tmp})
        ELSEIF("${_toggle}" STREQUAL "optimized")
          LIST(APPEND D2K_LIBRARIES_RELEASE ${_tmp})
        ENDIF()
      ENDIF()
    ENDFOREACH()
  ENDIF()

  FOREACH(_var ${D2K_LIST_SUFFIXES})
    IF(NOT "${_var}" STREQUAL "LIBRARIES" AND DEFINED ${_feature}_${_var})
      LIST(APPEND D2K_${_var} ${${_feature}_${_var}})
    ENDIF()
  ENDFOREACH()

  FOREACH(_var ${D2K_STRING_SUFFIXES})
    IF(DEFINED ${_feature}_${_var})
      ADD_FLAGS(D2K_${_var} "${${_feature}_${_var}}")
    ENDIF()
  ENDFOREACH()

ENDMACRO()
