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
  ParsedFiniteElement<1,1> fe_builder11("ParsedFiniteElement<1,1>", "FE_Q(2)");
  ParsedFiniteElement<2,2> fe_builder22("ParsedFiniteElement<2,2>",
					"FESystem[FE_Q(2)^d]",
					"u,u");
  ParsedFiniteElement<2,3> fe_builder23("ParsedFiniteElement<2,3>",
                                        "FESystem[FE_Q(2)^d-FE_DGP(1)]",
					"u,u,p");
  ParsedFiniteElement<3,3> fe_builder33("ParsedFiniteElement<2,3>", "FE_DGQ(2)");

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  FiniteElement<1,1> *fe11 = fe_builder11();
  FiniteElement<2,2> *fe22 = fe_builder22();
  FiniteElement<2,3> *fe23 = fe_builder23();
  FiniteElement<3,3> *fe33 = fe_builder33();

  deallog << "Generated fe11: " <<  fe11->get_name() << std::endl
          << "Generated fe22: " <<  fe22->get_name() << std::endl
          << "Generated fe23: " <<  fe23->get_name() << std::endl
          << "Generated fe33: " <<  fe33->get_name() << std::endl;
}
