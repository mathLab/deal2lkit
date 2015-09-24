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
# A small wrapper around ADD_LIBRARY that will define a target for each
# build type specified in D2K_BUILD_TYPES
#
# It is assumed that the desired compilation configuration is set via
#   D2K_LINKER_FLAGS_${build}
#   D2K_CXX_FLAGS_${build}
#   D2K_DEFINITIONS_${build}
#
# as well as the global (for all build types)
#   D2K_LINKER_FLAGS
#   D2K_CXX_FLAGS
#   D2K_DEFINITIONS
#

MACRO(D2K_ADD_LIBRARY _library)

  FOREACH(_build ${D2K_BUILD_TYPES})
    STRING(TOLOWER ${_build} _build_lowercase)

    ADD_LIBRARY(${_library}.${_build_lowercase}
      ${ARGN}
      )

    SET_TARGET_PROPERTIES(${_library}.${_build_lowercase} PROPERTIES
      LINK_FLAGS "${D2K_LINKER_FLAGS} ${D2K_LINKER_FLAGS_${_build}}"
      COMPILE_DEFINITIONS "${D2K_DEFINITIONS};${D2K_DEFINITIONS_${_build}}"
      COMPILE_FLAGS "${D2K_CXX_FLAGS} ${D2K_CXX_FLAGS_${_build}}"
      LINKER_LANGUAGE "CXX"
      )

    SET_PROPERTY(GLOBAL APPEND PROPERTY D2K_OBJECTS_${_build}
      "$<TARGET_OBJECTS:${_library}.${_build_lowercase}>"
      )

		#STRING(TOUPPER "${_build_type}" _BUILD_TYPE)
		#DEAL_II_SETUP_TARGET(${_library}.${_build_lowercase} ${_BUILD_TYPE})
		DEAL_II_SETUP_TARGET(${_library}.${_build_lowercase} ${_build})
  ENDFOREACH()

ENDMACRO()
