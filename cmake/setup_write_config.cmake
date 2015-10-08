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
## Copyright (C) 2014 - 2015 by the deal.II authors
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



########################################################################
#                                                                      #
#                Query for git repository information:                 #
#                                                                      #
########################################################################

D2K_QUERY_GIT_INFORMATION()

FILE(WRITE ${CMAKE_BINARY_DIR}/revision.log
  "###
  #
  #  Git information:
  #        Branch:   ${D2K_GIT_BRANCH}
  #        Revision: ${D2K_GIT_REVISION}
  #
  ###"
  )


########################################################################
#                                                                      #
#              Write a nice configuration summary to file:             #
#                                                                      #
########################################################################

SET(_log_detailed "${CMAKE_BINARY_DIR}/detailed.log")
SET(_log_summary  "${CMAKE_BINARY_DIR}/summary.log")
FILE(REMOVE ${_log_detailed} ${_log_summary})

MACRO(_both)
  # Write to both log files:
  FILE(APPEND ${_log_detailed} "${ARGN}")
  FILE(APPEND ${_log_summary} "${ARGN}")
ENDMACRO()
MACRO(_detailed)
  # Only write to detailed.log:
  FILE(APPEND ${_log_detailed} "${ARGN}")
ENDMACRO()
MACRO(_summary)
  # Only write to summary.log:
  FILE(APPEND ${_log_summary} "${ARGN}")
ENDMACRO()

_both(
"###
#
#  ${D2K_PACKAGE_NAME} configuration:
#
#        CMAKE_BUILD_TYPE:       ${CMAKE_BUILD_TYPE}
#        CMAKE_INSTALL_PREFIX:   ${CMAKE_INSTALL_PREFIX}
#        CMAKE_SOURCE_DIR:       ${CMAKE_SOURCE_DIR}
"
)
IF("${D2K_GIT_SHORTREV}" STREQUAL "")
  _both("#                                (version ${D2K_PACKAGE_VERSION})\n")
ELSE()
  _both("#                                (version ${D2K_PACKAGE_VERSION}, shortrev ${D2K_GIT_SHORTREV})\n")
ENDIF()
_both(
"#        CMAKE_BINARY_DIR:       ${CMAKE_BINARY_DIR}
#        CMAKE_CXX_COMPILER:     ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}
#                                ${CMAKE_CXX_COMPILER}
"
)

IF(CMAKE_C_COMPILER_WORKS)
  _detailed("#        CMAKE_C_COMPILER:       ${CMAKE_C_COMPILER}\n")
ENDIF()
IF(CMAKE_Fortran_COMPILER_WORKS)
  _detailed("#        CMAKE_Fortran_COMPILER: ${CMAKE_Fortran_COMPILER}\n")
ENDIF()
_detailed("#        CMAKE_GENERATOR:        ${CMAKE_GENERATOR}\n")

_both("#\n")

_detailed(
"#  Base configuration (prior to feature configuration):
#        D2K_CXX_FLAGS:            ${BASE_CXX_FLAGS}
"
)
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        D2K_CXX_FLAGS_RELEASE:    ${BASE_CXX_FLAGS_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        D2K_CXX_FLAGS_DEBUG:      ${BASE_CXX_FLAGS_DEBUG}\n")
ENDIF()

_detailed("#        D2K_LINKER_FLAGS:         ${BASE_LINKER_FLAGS}\n")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        D2K_LINKER_FLAGS_RELEASE: ${BASE_LINKER_FLAGS_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        D2K_LINKER_FLAGS_DEBUG:   ${BASE_LINKER_FLAGS_DEBUG}\n")
ENDIF()

_detailed("#        D2K_DEFINITIONS:          ${BASE_DEFINITIONS}\n")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        D2K_DEFINITIONS_RELEASE:  ${BASE_DEFINITIONS_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        D2K_DEFINITIONS_DEBUG:    ${BASE_DEFINITIONS_DEBUG}\n")
ENDIF()

_detailed("#        D2K_USER_DEFINITIONS:     ${BASE_DEFINITIONS}\n")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        D2K_USER_DEFINITIONS_REL: ${BASE_DEFINITIONS_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        D2K_USER_DEFINITIONS_DEB: ${BASE_DEFINITIONS_DEBUG}\n")
ENDIF()

_detailed("#        D2K_INCLUDE_DIRS          ${BASE_INCLUDE_DIRS}\n")
_detailed("#        D2K_USER_INCLUDE_DIRS:    ${BASE_USER_INCLUDE_DIRS}\n")
_detailed("#        D2K_BUNDLED_INCLUDE_DIRS: ${BASE_BUNDLED_INCLUDE_DIRS}\n")

_detailed("#        D2K_LIBRARIES:            ${BASE_LIBRARIES}\n")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        D2K_LIBRARIES_RELEASE:    ${BASE_LIBRARIES_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        D2K_LIBRARIES_DEBUG:      ${BASE_LIBRARIES_DEBUG}\n")
ENDIF()

_detailed("#\n")

_both("#  Configured Features:\n#\n")
GET_CMAKE_PROPERTY(_variables VARIABLES)
FOREACH(_var ${_variables})
  IF(_var MATCHES "D2K_WITH")
    LIST(APPEND _features "${_var}")
  ELSEIF(_var MATCHES "D2K_COMPONENT")
    LIST(APPEND _components "${_var}")
  ENDIF()
  IF(_var MATCHES "D2K_ENABLE_TESTING")
    LIST(APPEND _components "${_var}")
  ENDIF()

ENDFOREACH()

FOREACH(_var ${_features})
  IF(${${_var}})

    #
    # The feature is enabled:
    #
    STRING(REGEX REPLACE "^D2K_WITH_" "" _feature ${_var})
    IF(FEATURE_${_feature}_EXTERNAL_CONFIGURED)
      _both("#        ${_var} set up with external dependencies\n")
    ELSEIF(FEATURE_${_feature}_BUNDLED_CONFIGURED)
      IF(D2K_FORCE_BUNDLED_${_feature})
        _both("#        ${_var} set up with bundled packages (forced)\n")
      ELSE()
        _both("#        ${_var} set up with bundled packages\n")
      ENDIF()
    ELSE()
      _both("#        ${_var} = ${${_var}}\n")
    ENDIF()

    #
    # Print out version number:
    #
    IF(DEFINED ${_feature}_VERSION)
      _detailed("#            ${_feature}_VERSION = ${${_feature}_VERSION}\n")
    ENDIF()

    #
    # Print out ${_feature}_DIR:
    #
    IF(NOT "${${_feature}_DIR}" STREQUAL "")
      _detailed("#            ${_feature}_DIR = ${${_feature}_DIR}\n")
    ENDIF()

    #
    # Print the feature configuration:
    #
    FOREACH(_var2
        C_COMPILER CXX_COMPILER Fortran_COMPILER
        ${D2K_STRING_SUFFIXES} ${D2K_LIST_SUFFIXES}
        )
      IF(DEFINED ${_feature}_${_var2})
        _detailed("#            ${_feature}_${_var2} = ${${_feature}_${_var2}}\n")
      ENDIF()
    ENDFOREACH()
  ELSE()
    # FEATURE is disabled
    _both("#      ( ${_var} = ${${_var}} )\n")
  ENDIF()
ENDFOREACH()

_both(
"#\n#  Component configuration:\n#\n"
)
FOREACH(_var ${_components})
  IF(_var MATCHES "D2K_COMPONENT")
    IF(${${_var}})
      STRING(REPLACE "D2K_COMPONENT_" "" _component ${_var})
      LIST(APPEND _components ${_component})
      _both("#        ${_var} = ${${_var}}  \n")
    ELSE()
      _both("#      ( ${_var} = ${${_var}} )\n")
    ENDIF()
  ENDIF()
  IF(_var MATCHES "D2K_ENABLE_TESTING")
    IF(${${_var}})
      _both("#        ${_var} = ${${_var}}  \n")
    ELSE()
      _both("#      ( ${_var} = ${${_var}} )\n")
    ENDIF()
  ENDIF()
ENDFOREACH()

_summary(
"#\n#  Detailed information (compiler flags, feature configuration) can be found in detailed.log
#\n#  Run  $ "
)
IF(CMAKE_GENERATOR MATCHES "Ninja")
  _summary("ninja ")
ELSE()
  _summary("make ")
ENDIF()
_summary("info  to print a help message with a list of top level targets\n")

_both("#\n###")
