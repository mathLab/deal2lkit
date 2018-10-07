//-----------------------------------------------------------
//
//    Copyright (C) 2015 by the deal2lkit authors
//
//    This file is part of the deal2lkit library.
//
//    The deal2lkit library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal2lkit distribution.
//
//-----------------------------------------------------------


#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_finite_element.h>
#include <deal2lkit/utilities.h>

#include "../tests.h"


using namespace deal2lkit;


int
main()
{
  initlog();
  ParsedFiniteElement<1, 1> fe_builder11("ParsedFiniteElement<1,1>", "FE_Q(2)");
  ParsedFiniteElement<2, 2> fe_builder22("ParsedFiniteElement<2,2>",
                                         "FESystem[FE_Q(2)^d]",
                                         "u,u");
  ParsedFiniteElement<2, 3> fe_builder23("ParsedFiniteElement<2,3>",
                                         "FESystem[FE_Q(2)^d-FE_DGP(1)]",
                                         "u,u,p");
  ParsedFiniteElement<3, 3> fe_builder33("ParsedFiniteElement<2,3>",
                                         "FE_DGQ(2)");

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  FiniteElement<1, 1> *fe11 = fe_builder11();
  FiniteElement<2, 2> *fe22 = fe_builder22();
  FiniteElement<2, 3> *fe23 = fe_builder23();
  FiniteElement<3, 3> *fe33 = fe_builder33();

  deallog << "Generated fe11: " << fe11->get_name() << std::endl
          << "Generated fe22: " << fe22->get_name() << std::endl
          << "Generated fe23: " << fe23->get_name() << std::endl
          << "Generated fe33: " << fe33->get_name() << std::endl;
}
