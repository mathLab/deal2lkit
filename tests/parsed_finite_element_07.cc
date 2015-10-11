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


// Test if we correctly fail when the number of blocks is inconsistent

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
