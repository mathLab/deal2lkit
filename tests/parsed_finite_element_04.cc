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
#include "utilities.h"
#include "parameter_acceptor.h"
#include "parsed_finite_element.h"


int main ()
{
  initlog();
  ParsedFiniteElement<2> fe_builder("FE", "FESystem[FE_Q(2)^d-FE_DGP(1)]", "u,u,p", 1);

  ParameterAcceptor::initialize("test_in.prm", "test_out.prm");
  // Should fail because fe_name has 3 components
  return 0;
}
