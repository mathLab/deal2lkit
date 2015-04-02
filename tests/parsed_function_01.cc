//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

#include "tests.h"
#include "parsed_function.h"

int main ()
{
  initlog();

  ParsedFunction<2> pf;

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);


}
