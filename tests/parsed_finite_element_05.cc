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
#include <deal2lkit/parsed_finite_element.h>

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
