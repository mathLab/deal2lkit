##
# CMake script for large projects, with one executable
# and a library
##

# Set the name of the project and target library:
SET(TARGET_LIB dealii-sak)

# The source directory contains all implementations of the various classes.
# No main functions are allowed in these files.
file(GLOB _library_source_files source/*cc)

# All files in apps are assumed to contain a main function.
file(GLOB _applications_source_files apps/*cc)

INCLUDE_DIRECTORIES(include)

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.3 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET_LIB})

IF(DEAL_II_BUILD_TYPE MATCHES "Release")
    ADD_LIBRARY(${TARGET_LIB} SHARED ${_library_source_files})
    FOREACH(_main ${_applications_source_files})
	GET_FILENAME_COMPONENT(_target ${_main} NAME_WE)
	ADD_EXECUTABLE(${_target} ${_main})
	TARGET_LINK_LIBRARIES(${_target} ${TARGET_LIB})
	DEAL_II_SETUP_TARGET(${_target} RELEASE)
    ENDFOREACH()
    SET(TEST_LIBRARIES_RELEASE ${TARGET_LIB})
    DEAL_II_SETUP_TARGET(${TARGET_LIB} RELEASE)
ENDIF()

IF(DEAL_II_BUILD_TYPE MATCHES "Debug")
    ADD_LIBRARY(${TARGET_LIB}.g SHARED ${_library_source_files})   
    FOREACH(_main ${_applications_source_files})
	GET_FILENAME_COMPONENT(_target ${_main} NAME_WE)
	ADD_EXECUTABLE(${_target}.g ${_main})
	TARGET_LINK_LIBRARIES(${_target}.g ${TARGET_LIB}.g)
	DEAL_II_SETUP_TARGET(${_target}.g DEBUG)
    ENDFOREACH()
    SET(TEST_LIBRARIES_DEBUG ${TARGET_LIB}.g)
    DEAL_II_SETUP_TARGET(${TARGET_LIB}.g DEBUG)
ENDIF()

ADD_SUBDIRECTORY(tests)
ENABLE_TESTING()
