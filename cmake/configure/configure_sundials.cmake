## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2015 by the deal2lkit authors
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
# Configuration for the sundials library:
#

MACRO(FEATURE_SUNDIALS_FIND_EXTERNAL var)

	FIND_PACKAGE(SUNDIALS)

		IF(SUNDIALS_FOUND)

		SET(${var} TRUE)

    #
    # Set SUNDIALS_DIR to something meaningful if empty
    #
    IF("${SUNDIALS_DIR}" STREQUAL "")
      SET(SUNDIALS_DIR "<system location>")
    ENDIF()

	ENDIF()
ENDMACRO()

CONFIGURE_FEATURE(SUNDIALS)
