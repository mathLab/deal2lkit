#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

# If the result of a find_* instruction is -NOTFOUND, then 
# print an error message and exit. If we pass two arguments,
# then set the second variable name to TRUE or FALSE. In this case
# we do not exit if the result is -NOTFOUND.
#
MACRO(CHECK_FOUND _variable)
  IF("${${_variable}}" STREQUAL "${_variable}-NOTFOUND")
    IF("${ARGC}" STREQUAL "1")
	MESSAGE(FATAL_ERROR "Could not find ${_variable}!")
    ELSE()
	MESSAGE("Could not find ${_variable}. Setting ${ARGV1} to FALSE.")
	SET(${ARGV1} FALSE)
    ENDIF()
  ELSE()
    MESSAGE("Found ${_variable}: ${${_variable}}")
    IF("${ARGC}" STREQUAL "2")
	SET(${ARGV1} TRUE)
	MESSAGE("Setting ${ARGV1}=${${ARGV1}}.")
    ENDIF()
  ENDIF()
ENDMACRO()

