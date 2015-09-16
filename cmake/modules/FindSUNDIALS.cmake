## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2014 by the deal2lkit authors
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

#
# Try to find the SUNDIALS libraries
#
# This module exports
#
#   NETCDF_LIBRARIES
#   NETCDF_INCLUDE_DIRS
#

SET(SUNDIALS_DIR "" CACHE PATH "An optional hint to a SUNDIALS_DIR installation")
SET_IF_EMPTY(SUNDIALS_DIR "$ENV{SUNDIALS_DIR}")

D2K_FIND_LIBRARY(SUNDIALS_LIB_IDA NAMES sundials_ida
	HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib
  )

D2K_FIND_LIBRARY(SUNDIALS_LIB_SER NAMES sundials_nvecserial
	HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib
  )

D2K_FIND_LIBRARY(SUNDIALS_LIB_PAR NAMES sundials_nvecparallel
	HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib
  )

SET(SUN_INC "${SUNDIALS_DIR}/include")

INCLUDE_DIRECTORIES(${SUN_INC})

SET(SUNDIALS_LIBS "sundials_ida;sundials_nvecserial;sundials_nvecparallel")
FOREACH(_BUILD_TYPE ${D2K_BUILD_TYPES})
	SET(_lib ${D2K_BASE_NAME}${D2K_${_BUILD_TYPE}_SUFFIX})
	FOREACH(_slib ${SUNDIALS_LIBS})
		FIND_LIBRARY(SUNDIALS_${_slib} ${_slib} 
			HINTS ${SUNDIALS_DIR} 
			PATH_SUFFIXES lib)
		IF(NOT "${SUNDIALS_${_slib}}" STREQUAL "${_slib}-NOTFOUND")
			TARGET_LINK_LIBRARIES(${_lib} ${SUNDIALS_${_slib}})
			MESSAGE("-- Library ${_slib} found: ${SUNDIALS_${_slib}}")
		ELSE()
			MESSAGE(WARNING "-- Library ${_slib} not found: ${SUNDIALS_${_slib}}")
		ENDIF()
	ENDFOREACH()
ENDFOREACH()



D2K_PACKAGE_HANDLE(SUNDIALS
  LIBRARIES REQUIRED 
	SUNDIALS_LIB_IDA
	SUNDIALS_LIB_SER
	SUNDIALS_LIB_PAR
	INCLUDE_DIRS 
	REQUIRED SUN_INC
	USER_INCLUDE_DIRS
	REQUIRED SUN_INC
	 CLEAR 
	SUNDIALS_LIB_IDA
	SUNDIALS_LIB_SER
	SUNDIALS_LIB_PAR
	SUN_INC
  )

