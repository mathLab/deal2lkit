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
## Copyright (C) 2013 - 2015 by the deal.II authors
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
# Add convenience targets that build and install only a specific component:
#
#   library
#   documentation
#   examples
#


IF("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local")
  #
  # In case that CMAKE_INSTALL_PREFIX wasn't set, we assume that the user
  # doesn't actually want to install but just use deal2lkit in the build
  # directory. In this case, do not add the "install" phase to the
  # convenience targets.
  #
  MACRO(_add_custom_target _name)
    ADD_CUSTOM_TARGET(${_name})
  ENDMACRO()

  # Print precise informations about the convenience targets:
  SET(_description_string "build")
ELSE()
  MACRO(_add_custom_target _name)
    ADD_CUSTOM_TARGET(${_name}
      COMMAND ${CMAKE_COMMAND}
      -DCOMPONENT="${_name}" -P cmake_install.cmake
      COMMENT "Build and install component \"library\"."
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      )
  ENDMACRO()

  # Print precise informations about the convenience targets:
  SET(_description_string "build and install")
ENDIF()

# The library can always be compiled and/or installed unconditionally ;-)
_add_custom_target(library)

FOREACH(_component documentation examples)
  STRING(TOUPPER "${_component}" _component_uppercase)
  IF(D2K_COMPONENT_${_component_uppercase})
    _add_custom_target(${_component})
  ELSE()
    STRING(TOUPPER ${_component} _componentuppercase)
    ADD_CUSTOM_TARGET(${_component}
      COMMAND
      ${CMAKE_COMMAND} -E echo ''
      && ${CMAKE_COMMAND} -E echo ''
      && ${CMAKE_COMMAND} -E echo '***************************************************************************'
      && ${CMAKE_COMMAND} -E echo "**  Error: Could not ${_description_string} disabled component \"${_component}\"."
      && ${CMAKE_COMMAND} -E echo "**  Please reconfigure with -DD2K_COMPONENT_${_componentuppercase}=ON"
      && ${CMAKE_COMMAND} -E echo '***************************************************************************'
      && ${CMAKE_COMMAND} -E echo ''
      && ${CMAKE_COMMAND} -E echo ''
      && false
      )
  ENDIF()
ENDFOREACH()

IF(NOT D2K_COMPONENT_PACKAGE)
  ADD_CUSTOM_TARGET(package
    COMMAND
    ${CMAKE_COMMAND} -E echo ''
    && ${CMAKE_COMMAND} -E echo ''
    && ${CMAKE_COMMAND} -E echo '***************************************************************************'
    && ${CMAKE_COMMAND} -E echo "**  Error: Could not generate binary package. The component is disabled."
    && ${CMAKE_COMMAND} -E echo "**  Please reconfigure with -DD2K_COMPONENT_PACKAGE=ON"
    && ${CMAKE_COMMAND} -E echo '***************************************************************************'
    && ${CMAKE_COMMAND} -E echo ''
    && ${CMAKE_COMMAND} -E echo ''
    && false
    )
ENDIF()


ADD_CUSTOM_TARGET(indent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./scripts/indent
  )
