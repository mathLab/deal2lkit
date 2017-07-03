## ---------------------------------------------------------------------
##
## Copyright (C) 2015-2016 by the deal2lkit authors
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
#   SUNDIALS_LIBRARIES
#   SUNDIALS_INCLUDE_DIRS
#

SET(SUNDIALS_DIR "" CACHE PATH "An optional hint to a SUNDIALS_DIR installation")
SET_IF_EMPTY(SUNDIALS_DIR "$ENV{SUNDIALS_DIR}")

D2K_FIND_LIBRARY(SUNDIALS_LIB_IDA NAMES sundials_ida
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib
  )

D2K_FIND_LIBRARY(SUNDIALS_LIB_KINSOL NAMES sundials_kinsol
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib
  )

 D2K_FIND_LIBRARY(SUNDIALS_LIB_ARKODE NAMES sundials_arkode
   HINTS ${SUNDIALS_DIR}
   PATH_SUFFIXES lib
   )

D2K_FIND_LIBRARY(SUNDIALS_LIB_SER NAMES sundials_nvecserial
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib
  )

SET(SUN_INC "${SUNDIALS_DIR}/include")

IF(DEAL_II_WITH_MPI)
  D2K_FIND_LIBRARY(SUNDIALS_LIB_PAR NAMES sundials_nvecparallel
    HINTS ${SUNDIALS_DIR}
    PATH_SUFFIXES lib
    )

  D2K_PACKAGE_HANDLE(SUNDIALS
    LIBRARIES REQUIRED
    SUNDIALS_LIB_IDA
    SUNDIALS_LIB_KINSOL
    SUNDIALS_LIB_ARKODE
    SUNDIALS_LIB_SER
    SUNDIALS_LIB_PAR
    INCLUDE_DIRS
    REQUIRED SUN_INC
    USER_INCLUDE_DIRS
    REQUIRED SUN_INC
    CLEAR
    SUNDIALS_LIB_IDA
    SUNDIALS_LIB_KINSOL
    SUNDIALS_LIB_ARKODE
    SUNDIALS_LIB_SER
    SUNDIALS_LIB_PAR
    SUN_INC
    )

ELSE()

  D2K_PACKAGE_HANDLE(SUNDIALS
    LIBRARIES REQUIRED
    SUNDIALS_LIB_IDA
    SUNDIALS_LIB_KINSOL
    SUNDIALS_LIB_ARKODE
    SUNDIALS_LIB_SER
    INCLUDE_DIRS
    REQUIRED SUN_INC
    USER_INCLUDE_DIRS
    REQUIRED SUN_INC
    CLEAR
    SUNDIALS_LIB_IDA
    SUNDIALS_LIB_KINSOL
    SUNDIALS_LIB_ARKODE
    SUNDIALS_LIB_SER
    SUN_INC
    )

ENDIF()
