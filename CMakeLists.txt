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

DEAL_II_INITIALIZE_CACHED_VARIABLES()
DEAL_II_QUERY_GIT_INFORMATION()
PROJECT(${TARGET_LIB})

# We generate one library for each type of deal.II library
# we found. If you compiled deal.II with both Release and Debug
# mode, this will generate both Release and Debug programs for you
# The debug library and program are postfixed with ".g"
SET(_d2_build_types "Release;Debug")
SET(Release_postfix "")
SET(Debug_postfix ".g")


FOREACH(_build_type ${_d2_build_types})
    # Only build this type, if deal.II was compiled with it.
    IF(DEAL_II_BUILD_TYPE MATCHES "${_build_type}")
	MESSAGE("-- Found ${_build_type} version of deal.II.")
	SET(_postfix "${${_build_type}_postfix}")

	SET(_lib "${TARGET_LIB}${_postfix}")
	MESSAGE("-- Configuring library ${_lib}")

	STRING(TOUPPER "${_build_type}" _BUILD_TYPE)

	# Add all files, except the main file, to the library
	ADD_LIBRARY(${_lib} SHARED ${_library_source_files})
        INSTALL(TARGETS ${_lib} LIBRARY DESTINATION lib)

	SET_PROPERTY(TARGET ${_lib} APPEND PROPERTY COMPILE_DEFINITIONS
	  DEAL_II_SAK_GIT_BRANCH="${GIT_BRANCH}${_postfix}"
	  DEAL_II_SAK_GIT_REVISION="${GIT_REVISION}${_postfix}"
	  DEAL_II_SAK_GIT_SHORTREV="${GIT_SHORTREV}${_postfix}"
	  DEAL_II_VERSION="${DEAL_II_VERSION}${_postfix}"
	  DEAL_II_GIT_BRANCH="${DEAL_II_GIT_BRANCH}${_postfix}"
	  DEAL_II_GIT_REVISION="${DEAL_II_GIT_REVISION}${_postfix}"
	  DEAL_II_GIT_SHORTREV="${DEAL_II_GIT_SHORTREV}${_postfix}"
	)
    
	# Compile an executable only if required by the user

    	FOREACH(_main ${_applications_source_files})
		# Find out the the name of the executable
    		GET_FILENAME_COMPONENT(_target_a ${_main} NAME_WE)
		SET(_target "${_target_a}${${_build_type}_postfix}")
		MESSAGE("-- Configuring executable ${_target}")

    		ADD_EXECUTABLE(${_target} ${_main})
    		INSTALL(TARGETS ${_target} RUNTIME DESTINATION bin)
    		TARGET_LINK_LIBRARIES(${_target} ${_lib})
		DEAL_II_SETUP_TARGET(${_target} ${_BUILD_TYPE})
    	ENDFOREACH()

	SET(TEST_LIBRARIES_${_BUILD_TYPE} ${_lib})
	DEAL_II_SETUP_TARGET(${_lib} ${_BUILD_TYPE})
    ENDIF()
ENDFOREACH()

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

