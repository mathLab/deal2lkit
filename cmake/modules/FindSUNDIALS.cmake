## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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
