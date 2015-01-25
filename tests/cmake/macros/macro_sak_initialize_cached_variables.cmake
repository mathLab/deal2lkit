#
# This file implements the SAK_INITIALIZE_VARIABLES macro, which is
# part of the SAK library.
#
# Usage:
#       SAK_INITIALIZE_CACHED_VARIABLES()
#
# This macro has to be called before PROJECT()!
#

MACRO(SAK_INITIALIZE_CACHED_VARIABLES)

  #
  # Set build type according to available libraries
  #
  IF(SAK_BUILD_TYPE MATCHES "Debug")
    SET(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug, Release"
      )
  ELSE()
    SET(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug, Release"
      )
  ENDIF()

  #
  # Bail out if build type is unknown...
  #
  IF( NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release" AND
      NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
    MESSAGE(FATAL_ERROR
      "\nCMAKE_BUILD_TYPE does neither match Release nor Debug!\n\n"
      )
  ENDIF()
  #
  # ... or unsupported
  #
  IF(NOT SAK_BUILD_TYPE MATCHES "${CMAKE_BUILD_TYPE}")
    MESSAGE(FATAL_ERROR "\n"
      "CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\" unsupported by current installation!\n"
      "SAK was build with \"${SAK_BUILD_TYPE}\" only build type.\n\n"
      )
  ENDIF()

  SET(CMAKE_CXX_COMPILER ${DEAL_II_CXX_COMPILER} CACHE STRING
    "CXX Compiler.")

  SET(CMAKE_CXX_FLAGS ${DEAL_II_CXX_FLAGS} CACHE STRING
    "Flags used by the compiler during all build types."
    )

  SET(CMAKE_CXX_FLAGS_DEBUG ${DEAL_II_CXX_FLAGS_DEBUG} CACHE STRING
    "Flags used by the compiler during debug builds."
    )

  SET(CMAKE_CXX_FLAGS_RELEASE ${DEAL_II_CXX_FLAGS_RELEASE} CACHE STRING
    "Flags used by the compiler during release builds."
    )

  MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)

ENDMACRO()

