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
    INSTALL(TARGETS ${TARGET_LIB} LIBRARY DESTINATION lib)
    FOREACH(_main ${_applications_source_files})
	GET_FILENAME_COMPONENT(_target ${_main} NAME_WE)
	ADD_EXECUTABLE(${_target} ${_main})
	INSTALL(TARGETS ${_target} RUNTIME DESTINATION bin)
	TARGET_LINK_LIBRARIES(${_target} ${TARGET_LIB})
	DEAL_II_SETUP_TARGET(${_target} RELEASE)
    ENDFOREACH()
    SET(TEST_LIBRARIES_RELEASE ${TARGET_LIB})
    DEAL_II_SETUP_TARGET(${TARGET_LIB} RELEASE)
ENDIF()

IF(DEAL_II_BUILD_TYPE MATCHES "Debug")
    ADD_LIBRARY(${TARGET_LIB}.g SHARED ${_library_source_files})   
    INSTALL(TARGETS ${TARGET_LIB}.g LIBRARY DESTINATION lib)
    FOREACH(_main ${_applications_source_files})
	GET_FILENAME_COMPONENT(_target ${_main} NAME_WE)
	ADD_EXECUTABLE(${_target}.g ${_main})
	INSTALL(TARGETS ${_target}.g RUNTIME DESTINATION bin)
	TARGET_LINK_LIBRARIES(${_target}.g ${TARGET_LIB}.g)
	DEAL_II_SETUP_TARGET(${_target}.g DEBUG)
    ENDFOREACH()
    SET(TEST_LIBRARIES_DEBUG ${TARGET_LIB}.g)
    DEAL_II_SETUP_TARGET(${TARGET_LIB}.g DEBUG)
ENDIF()

ADD_SUBDIRECTORY(tests)
ENABLE_TESTING()

# Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include
        FILES_MATCHING PATTERN "*.h")

