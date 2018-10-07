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
  ParsedFiniteElement<1, 1> fe_builder11;
  ParsedFiniteElement<1, 2> fe_builder12;
  ParsedFiniteElement<1, 3> fe_builder13;
  ParsedFiniteElement<2, 2> fe_builder22;
  ParsedFiniteElement<2, 3> fe_builder23;
  ParsedFiniteElement<3, 3> fe_builder33;

  ParameterHandler prm;
  ParameterAcceptor::declare_all_parameters(prm);
  ParameterAcceptor::parse_all_parameters(prm);
  prm.log_parameters(deallog);

  auto fe11 = fe_builder11();
  auto fe12 = fe_builder12();
  auto fe13 = fe_builder13();
  auto fe22 = fe_builder22();
  auto fe23 = fe_builder23();
  auto fe33 = fe_builder33();

  deallog << "Generated fe11: " << type(*fe11) << std::endl
          << "Generated fe12: " << type(*fe12) << std::endl
          << "Generated fe13: " << type(*fe13) << std::endl
          << "Generated fe22: " << type(*fe22) << std::endl
          << "Generated fe23: " << type(*fe23) << std::endl
          << "Generated fe33: " << type(*fe33) << std::endl;
}
