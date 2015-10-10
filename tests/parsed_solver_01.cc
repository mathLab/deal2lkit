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
#include <deal2lkit/parsed_solver.h>

#include <deal.II/lac/vector.h>


using namespace deal2lkit;

int main ()
{
  initlog();

  ParsedSolver<Vector<double> > solver("Solver", "cg", 100, 1e-6);
  ParameterAcceptor::initialize();

  ParameterAcceptor::prm.log_parameters(deallog);

  Vector<double> b(4);
  for (unsigned int i=0; i<4; ++i)
    {
      b(i) = i+1.;
    }
  Vector<double> x(4);

  solver.vmult(x,b);

  x -= b;

  deallog << "Norm: " << x.l2_norm() << std::endl;

}
