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


LIST(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/tests/cmake/
  ${CMAKE_SOURCE_DIR}/tests/cmake/modules/
  )

#
# Load all macros:
#
INCLUDE(setup_external_macros)
FILE(GLOB _macro_files "${CMAKE_SOURCE_DIR}/tests/cmake/macros/*.cmake")
FOREACH(_file ${_macro_files})
  INCLUDE(${_file})
ENDFOREACH()

#
# Setup the testsuite.
#
# We define toplevel targets:
#    setup_tests    - set up testsuite subprojects
#    regen_tests    - rerun configure stage in every testsuite subproject
#    clean_tests    - run the 'clean' target in every testsuite subproject
#    prune_tests    - remove all testsuite subprojects
#

SET_IF_EMPTY(MAKEOPTS $ENV{MAKEOPTS})

# Write minimalistic CTestTestfile.cmake files to CMAKE_BINARY_DIR and
# CMAKE_BINARY_DIR/tests:
#

MESSAGE(STATUS "Setting up testsuite")

FILE(WRITE ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "SUBDIRS(tests)")
FILE(WRITE ${CMAKE_BINARY_DIR}/tests/CTestTestfile.cmake "")

#
# A config line for tests:
#

SET(_options "")
LIST(APPEND _options -DSAK_SOURCE_DIR=${CMAKE_SOURCE_DIR})
LIST(APPEND _options -DSAK_BINARY_DIR=${CMAKE_BINARY_DIR})
LIST(APPEND _options -DSAK_TARGET=${TARGET})
LIST(APPEND _options -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
FOREACH(_var
  DIFF_DIR NUMDIFF_DIR TEST_DIFF TEST_PICKUP_REGEX TEST_TIME_LIMIT
  MPIEXEC MPIEXEC_NUMPROC_FLAG MPIEXEC_PREFLAGS MPIEXEC_POSTFLAGS
  )
  LIST(APPEND _options "-U${_var}")
  IF(DEFINED ${_var})
    LIST(APPEND _options "-D${_var}=${${_var}}")
  ENDIF()
ENDFOREACH()

SET(_categories
   template_category
  )

#
# Custom targets for the testsuite:
#

# Setup tests:
ADD_CUSTOM_TARGET(setup_tests)
ADD_DEPENDENCIES(setup_tests ${TARGET}lib)

# Remove all tests:
ADD_CUSTOM_TARGET(prune_tests)

# Regenerate tests (run "make rebuild_cache" in subprojects):
ADD_CUSTOM_TARGET(regen_tests)

# Regenerate tests (run "make clean" in subprojects):
ADD_CUSTOM_TARGET(clean_tests)

FOREACH(_category ${_categories})
  SET(_category_dir ${CMAKE_SOURCE_DIR}/tests/${_category})

  ADD_CUSTOM_TARGET(setup_tests_${_category}
    COMMAND ${CMAKE_COMMAND} -E make_directory
      ${CMAKE_BINARY_DIR}/tests/${_category}
    COMMAND cd ${CMAKE_BINARY_DIR}/tests/${_category} &&
      ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} ${_options} ${_category_dir}
      > /dev/null
    COMMENT "Processing tests/${_category}"
    )
  ADD_DEPENDENCIES(setup_tests setup_tests_${_category})

  ADD_CUSTOM_TARGET(prune_tests_${_category}
    COMMAND ${CMAKE_COMMAND} -E remove_directory
      ${CMAKE_BINARY_DIR}/tests/${_category}
    COMMENT "Processing tests/${_category}"
    )
  ADD_DEPENDENCIES(prune_tests prune_tests_${_category})

  ADD_CUSTOM_TARGET(regen_tests_${_category}
    COMMAND
      test ! -d ${CMAKE_BINARY_DIR}/tests/${_category} || ${CMAKE_COMMAND}
      --build ${CMAKE_BINARY_DIR}/tests/${_category} --target regenerate
    COMMENT "Processing tests/${_category}"
    )
  ADD_DEPENDENCIES(regen_tests regen_tests_${_category})

  ADD_CUSTOM_TARGET(clean_tests_${_category}
    COMMAND
      test ! -d ${CMAKE_BINARY_DIR}/tests/${_category} || ${CMAKE_COMMAND}
      --build ${CMAKE_BINARY_DIR}/tests/${_category} --target clean
    COMMENT "Processing tests/${_category}"
    )
  ADD_DEPENDENCIES(clean_tests clean_tests_${_category})

  FILE(APPEND ${CMAKE_BINARY_DIR}/tests/CTestTestfile.cmake
    "SUBDIRS(${_category})\n"
    )
ENDFOREACH()

MESSAGE(STATUS "Setting up testsuite - Done")

