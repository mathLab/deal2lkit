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

#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_finite_element.h>

int main ()
{
  initlog();
  ParsedFiniteElement<1,1> fe_builder11;
  ParsedFiniteElement<1,2> fe_builder12;
  ParsedFiniteElement<1,3> fe_builder13;
  ParsedFiniteElement<2,2> fe_builder22;
  ParsedFiniteElement<2,3> fe_builder23;
  ParsedFiniteElement<3,3> fe_builder33;

  ParameterHandler prm;
  ParameterAcceptor::declare_all_parameters(prm);
  ParameterAcceptor::parse_all_parameters(prm);
  prm.log_parameters(deallog);

  FiniteElement<1,1> *fe11 = fe_builder11();
  FiniteElement<1,2> *fe12 = fe_builder12();
  FiniteElement<1,3> *fe13 = fe_builder13();
  FiniteElement<2,2> *fe22 = fe_builder22();
  FiniteElement<2,3> *fe23 = fe_builder23();
  FiniteElement<3,3> *fe33 = fe_builder33();

  deallog << "Generated fe11: " <<  type(*fe11) << std::endl
          << "Generated fe12: " <<  type(*fe12) << std::endl
          << "Generated fe13: " <<  type(*fe13) << std::endl
          << "Generated fe22: " <<  type(*fe22) << std::endl
          << "Generated fe23: " <<  type(*fe23) << std::endl
          << "Generated fe33: " <<  type(*fe33) << std::endl;
}
