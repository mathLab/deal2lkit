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

// Test if we correctly fail when the number of blocks is inconsistent

#include "tests.h"
#include "parsed_finite_element.h"

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

int main ()
{
  initlog();

  ParsedFiniteElement<2,2> fe("Finite Element",
                              "FESystem[FE_Q(2)^2-FE_DGP(1)^2]",
                              "u,u,w,w", 4, "1,0,2; 0,1,2");

  ParameterAcceptor::initialize();

  ParameterAcceptor::prm.log_parameters(deallog);

  Table<2, DoFTools::Coupling> c = fe.get_coupling();
  log(c);

}
