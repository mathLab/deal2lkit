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

MESSAGE(STATUS "Setting up library")

#
# Compile the deal.II library
#

INCLUDE_DIRECTORIES(
  ${CMAKE_BINARY_DIR}/include/
  ${CMAKE_SOURCE_DIR}/include/
  ${D2K_INCLUDE_DIRS}
  )

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_src
	utilities.cc
	parameter_acceptor.cc
	parsed_finite_element.cc
	parsed_grid_generator.cc
	parsed_data_out.cc
	error_handler.cc
	parsed_mapped_functions.cc
	dae_time_integrator.cc
	sundials_interface.cc
	)

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal2lkit/*.h
  )

D2K_ADD_LIBRARY(obj_base OBJECT ${_src} ${_header} ${_inst}
  ${CMAKE_BINARY_DIR}/include/deal2lkit/config.h
  )

FOREACH(build ${D2K_BUILD_TYPES})
  STRING(TOLOWER ${build} build_lowercase)

  #
  # Combine all ${build} OBJECT targets to a ${build} library:
  #

  GET_PROPERTY(_objects GLOBAL PROPERTY D2K_OBJECTS_${build})
  ADD_LIBRARY(${D2K_BASE_NAME}${D2K_${build}_SUFFIX}
    dummy.cc # Workaround for a bug in the Xcode generator
    ${_objects}
    )
  ADD_DEPENDENCIES(library ${D2K_BASE_NAME}${D2K_${build}_SUFFIX})

  SET_TARGET_PROPERTIES(${D2K_BASE_NAME}${D2K_${build}_SUFFIX}
    PROPERTIES
    VERSION "${D2K_PACKAGE_VERSION}"
    #
    # Sonaming: Well... we just use the version number.
    # No point to wrack one's brain over the question whether a new version of
    # a C++ library is still ABI backwards compatible :-]
    #
    SOVERSION "${D2K_PACKAGE_VERSION}"
    LINK_FLAGS "${D2K_LINKER_FLAGS} ${D2K_LINKER_FLAGS_${build}}"
    LINKER_LANGUAGE "CXX"
    COMPILE_DEFINITIONS "${D2K_DEFINITIONS};${D2K_DEFINITIONS_${build}}"
    COMPILE_FLAGS "${D2K_CXX_FLAGS} ${D2K_CXX_FLAGS_${build}}"
    INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${D2K_LIBRARY_RELDIR}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${D2K_LIBRARY_RELDIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${D2K_LIBRARY_RELDIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${D2K_EXECUTABLE_RELDIR}"
    )

  TARGET_LINK_LIBRARIES(${D2K_BASE_NAME}${D2K_${build}_SUFFIX}
    ${D2K_LIBRARIES_${build}}
    ${D2K_LIBRARIES}
    )

  FILE(MAKE_DIRECTORY
    ${CMAKE_BINARY_DIR}/${D2K_PROJECT_CONFIG_RELDIR}
    )
  EXPORT(TARGETS ${D2K_BASE_NAME}${D2K_${build}_SUFFIX}
    FILE ${CMAKE_BINARY_DIR}/${D2K_PROJECT_CONFIG_RELDIR}/${D2K_PROJECT_CONFIG_NAME}Targets.cmake
    APPEND
    )

  INSTALL(TARGETS ${D2K_BASE_NAME}${D2K_${build}_SUFFIX}
    COMPONENT library
    EXPORT ${D2K_PROJECT_CONFIG_NAME}Targets
    RUNTIME DESTINATION ${D2K_EXECUTABLE_RELDIR}
    LIBRARY DESTINATION ${D2K_LIBRARY_RELDIR}
    ARCHIVE DESTINATION ${D2K_LIBRARY_RELDIR}
    )
ENDFOREACH()

INSTALL(EXPORT ${D2K_PROJECT_CONFIG_NAME}Targets
  DESTINATION ${D2K_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )

MESSAGE(STATUS "Setting up library - Done")
