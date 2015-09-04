##
# Deal.II Swiss Army Knife
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

# Set standard flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-pedantic")
SET(CMAKE_BUILD_TYPE "${DEAL_II_BUILD_TYPE}" CACHE STRING "Build type. Either Debug, Release, or DebugRelease.")

DEAL_II_INITIALIZE_CACHED_VARIABLES()
DEAL_II_QUERY_GIT_INFORMATION()
PROJECT(${TARGET_LIB})


# We generate one library for each type of deal.II library
# we found. If you compiled deal.II with both Release and Debug
# mode, this will generate both Release and Debug programs for you
# The debug library and program are postfixed with ".g"

# Options about examples
OPTION(DEAL_II_SAK_ENABLE_TESTING "Enable testing." OFF)
OPTION(DEAL_II_SAK_COMPILE_EXAMPLES "Compile examples" OFF)
OPTION(DEAL_II_SAK_INSTALL_EXAMPLES "Install examples" ON)
OPTION(DEAL_II_SAK_USE_BOOST "Use boost libraries." ON)


IF(DEAL_II_SAK_USE_BOOST)
  FIND_PACKAGE(Boost COMPONENTS system filesystem)
  IF(Boost_FOUND)
    # add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})
    SET(Boost_USE_STATIC_LIBS OFF)
    SET(Boost_USE_MULTITHREADED ON)
    SET(Boost_USE_STATIC_RUNTIME OFF)
    ADD_DEFINITIONS(-DDEAL_II_SAK_WITH_BOOST)
    INCLUDE_DIRECTORIES (${BOOST_INCLUDE_DIRS})
  ENDIF()
ENDIF()

OPTION(DEAL_II_SAK_WITH_SUNDIALS "Use SUNDIALS libraries." OFF)
IF(DEAL_II_SAK_WITH_SUNDIALS)
  SET(SUNDIALS_DIR "$ENV{SUNDIALS_DIR}" CACHE PATH "Path to sundials")
  SET(SUNDIALS_LIBS "sundials_ida;sundials_nvecserial;sundials_nvecparallel")
  ADD_DEFINITIONS(-DDEAL_II_SAK_WITH_SUNDIALS)
  INCLUDE_DIRECTORIES(${SUNDIALS_DIR}/include)

  # the following option is useful to take advantage of the
  # macro of deal.II pick_up_test and allow us to
  # name "test.with_sundials=on.output" the output file of a test
  # in order run it only if sundials is installed.
  set (DEAL_II_WITH_SUNDIALS ON)
ELSE()
  set (DEAL_II_WITH_SUNDIALS OFF)
ENDIF()

# Loop over all build types.
SET(_d2_build_types "Release;Debug")
SET(Release_postfix "")
SET(Debug_postfix ".g")
FOREACH(_build_type ${_d2_build_types})
    # Only build this type, if deal.II was compiled with it.
    IF(DEAL_II_BUILD_TYPE MATCHES "${_build_type}" AND CMAKE_BUILD_TYPE MATCHES "${_build_type}")
	MESSAGE("-- Found ${_build_type} version of deal.II.")
	SET(_postfix "${${_build_type}_postfix}")

	SET(_lib "${TARGET_LIB}${_postfix}")
	MESSAGE("-- Configuring library ${_lib}")

	STRING(TOUPPER "${_build_type}" _BUILD_TYPE)

	# Add all files, except the main file, to the library
	ADD_LIBRARY(${_lib} SHARED ${_library_source_files})

  IF(DEAL_II_SAK_WITH_SUNDIALS)
    FOREACH(_slib ${SUNDIALS_LIBS})
      FIND_LIBRARY(SUNDIALS_${_slib} ${_slib} HINTS ${SUNDIALS_DIR} PATH_SUFFIXES lib)
      IF(NOT "${SUNDIALS_${_slib}}" STREQUAL "${_slib}-NOTFOUND")
        TARGET_LINK_LIBRARIES(${_lib} ${SUNDIALS_${_slib}})
        MESSAGE("-- Library ${_slib} found: ${SUNDIALS_${_slib}}")
      ELSE()
        MESSAGE(WARNING "-- Library ${_slib} not found: ${SUNDIALS_${_slib}}")
      ENDIF()
    ENDFOREACH()
  ENDIF()


  # Add boost dependencies if necessary
  IF(Boost_FOUND)
     TARGET_LINK_LIBRARIES(${_lib} ${Boost_LIBRARIES})
  ENDIF()
        INSTALL(TARGETS ${_lib} LIBRARY DESTINATION lib)

  SET_PROPERTY(TARGET ${_lib} APPEND PROPERTY COMPILE_DEFINITIONS
    DEAL_II_SAK_GIT_BRANCH="${GIT_BRANCH}${_postfix}"
    DEAL_II_SAK_GIT_REVISION="${GIT_REVISION}${_postfix}"
    DEAL_II_SAK_GIT_SHORTREV="${GIT_SHORTREV}${_postfix}"
  )

  # Compile an executable only if required by the user
  IF(DEAL_II_SAK_COMPILE_EXAMPLES)
         FOREACH(_main ${_applications_source_files})
    # Find out the the name of the executable
        GET_FILENAME_COMPONENT(_target_a ${_main} NAME_WE)
    SET(_target "${_target_a}${${_build_type}_postfix}")
    MESSAGE("-- Configuring executable ${_target}")

        ADD_EXECUTABLE(${_target} ${_main})
    IF(DEAL_II_SAK_INSTALL_EXAMPLES)
      INSTALL(TARGETS ${_target} RUNTIME DESTINATION bin)
    ENDIF()
        TARGET_LINK_LIBRARIES(${_target} ${_lib})
    DEAL_II_SETUP_TARGET(${_target} ${_BUILD_TYPE})
          ENDFOREACH()
  ENDIF()

  SET(TEST_LIBRARIES_${_BUILD_TYPE} ${_lib})
  DEAL_II_SETUP_TARGET(${_lib} ${_BUILD_TYPE})
    ENDIF()
ENDFOREACH()

IF(${DEAL_II_SAK_ENABLE_TESTING})
    ADD_SUBDIRECTORY(tests)
    ENABLE_TESTING()
ENDIF()

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
