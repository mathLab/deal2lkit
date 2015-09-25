## ---------------------------------------------------------------------
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
# A small macro to reset the CMAKE_REQUIRED_* variables to its default
# values
#
# Usage:
#     RESET_CMAKE_REQUIRED_FLAGS
#

MACRO(RESET_CMAKE_REQUIRED)
  SET(CMAKE_REQUIRED_FLAGS ${D2K_CXX_FLAGS_SAVED})
  SET(CMAKE_REQUIRED_INCLUDES)
  SET(CMAKE_REQUIRED_LIBRARIES ${D2K_LINKER_FLAGS_SAVED})
ENDMACRO()

