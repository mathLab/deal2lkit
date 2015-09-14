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
#include <deal2lkit/parsed_solver.h>

#include <deal.II/lac/vector.h>

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
