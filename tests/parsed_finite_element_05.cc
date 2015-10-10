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

#include "tests.h"
#include <deal2lkit/parsed_finite_element.h>


using namespace deal2lkit;

void log(const Table<2, DoFTools::Coupling> &c)
{
  for (unsigned int i=0; i<c.size(0); ++i)
    {
      deallog << "Row " << i << ": ";
      for (unsigned int j=0; j<c.size(1); ++j)
        deallog << (int) c[i][j] << ",";
      deallog << std::endl;
    }
}

void test(ParsedFiniteElement<2,2> &fe, const std::string &coupling,
          const std::string &prec_coupling="")
{
  ParameterAcceptor::prm.enter_subsection("Finite Element");
  ParameterAcceptor::prm.set("Block coupling", coupling);
  ParameterAcceptor::prm.set("Preconditioner block coupling", prec_coupling);
  ParameterAcceptor::prm.leave_subsection();
  ParameterAcceptor::parse_all_parameters();
  ParameterAcceptor::prm.log_parameters(deallog);

  Table<2, DoFTools::Coupling> c = fe.get_coupling();
  deallog << std::endl << std::endl
          << "Input string: " << coupling << std::endl;
  log(c);

  Table<2, DoFTools::Coupling> c2 = fe.get_preconditioner_coupling();
  deallog << std::endl << std::endl
          << "Input string: " << prec_coupling << std::endl;
  log(c2);
}


int main ()
{
  initlog();

  ParsedFiniteElement<2,2> fe("Finite Element",
                              "FESystem[FE_Q(2)^2-FE_DGP(1)]",
                              "u,u,p");

  ParameterAcceptor::declare_all_parameters();

  test(fe, "");
  test(fe, "1,1; 1,0");
  test(fe, "1,1,1; 1,1,1; 1,1,0");

  test(fe, "0,1; 1,0", "1,0; 0,1");

}
