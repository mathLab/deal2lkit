## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II authors
## Copyright (C) 2015 by the deal2lkit authors
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

IF(D2K_COMPONENT_EXAMPLES)
  MESSAGE(STATUS "Setting up examples")

  INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${D2K_EXAMPLES_RELDIR}
    COMPONENT examples
    FILES_MATCHING
    #
    # Exclude folder structures: doc, doxygen, CMakeFiles,...
    #
    PATTERN "CMakeFiles*" EXCLUDE
    PATTERN "doc*" EXCLUDE
    PATTERN "build*" EXCLUDE
    #
    # Glob Includes:
    #
    PATTERN "*.cc"
    PATTERN "*.h"
    PATTERN "*.prm"
    PATTERN "*.inp"
    PATTERN "*.output"
    #
    # All CMakeLists.txt except the one in examples/
    #
    PATTERN "assimp_importer/CMakeLists.txt"
    PATTERN "heat_equation/*/CMakeLists.txt"
    #
    # Special files:
    #
    #PATTERN "particular_file.special_format"
    )

  MESSAGE(STATUS "Setting up examples - Done")
ENDIF()
