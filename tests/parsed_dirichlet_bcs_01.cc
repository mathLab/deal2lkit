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
//test constructor
//
#include "tests.h"
#include "utilities.h"
#include "parameter_acceptor.h"
#include "parsed_dirichlet_bcs.h"

int main ()
{
  initlog();

  ParsedDirichletBCs<2,2,1> pf("Dirichlet","u","0=ALL","0=0");

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);


}
