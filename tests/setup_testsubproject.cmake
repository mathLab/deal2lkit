FIND_PACKAGE(deal.II 8.4 REQUIRED HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR})
FIND_PACKAGE(deal2lkit 1.0 REQUIRED HINTS ${CMAKE_BINARY_DIR}/../..)
DEAL_II_INITIALIZE_CACHED_VARIABLES()
D2K_INITIALIZE_CACHED_VARIABLES()
PROJECT(testsuite CXX)
SET(CMAKE_BUILD_TYPE ${D2K_BUILD_TYPE} CACHE STRING "" FORCE)

INCLUDE_DIRECTORIES(${D2K_INCLUDE_DIRS})

FOREACH(_var DIFF_DIR NUMDIFF_DIR TEST_PICKUP_REGEX TEST_TIME_LIMIT)
  SET_IF_EMPTY(${_var} "$ENV{${_var}}")
  SET(${_var} "${${_var}}" CACHE STRING "" FORCE)
ENDFOREACH()


# For each cofigured feature of deal2lkit we set the corresponding
# variable DEAL_II_WITH_featurename in order to exploit what is done
# for the deal.II library

GET_CMAKE_PROPERTY(_variables VARIABLES)
FOREACH(_var ${_variables})
  IF(_var MATCHES "D2K_WITH")
    LIST(APPEND _features "${_var}")
    ENDIF()
ENDFOREACH()

FOREACH(_var ${_features})
STRING(REGEX REPLACE "^D2K_WITH_" "" _feature ${_var})
SET(DEAL_II_WITH_${_feature} ${${_var}})

ENDFOREACH()

SET(TEST_LIBRARIES_DEBUG ${D2K_TARGET_DEBUG})
SET(TEST_LIBRARIES_RELEASE ${D2K_TARGET_RELEASE})

INCLUDE(${D2K_TARGET_CONFIG})

  
#
# A custom target that does absolutely nothing. It is used in the main
# project to trigger a "make rebuild_cache" if necessary.
#
ADD_CUSTOM_TARGET(regenerate)
