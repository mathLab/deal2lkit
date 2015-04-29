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
#include "parsed_function.h"

int main ()
{
  initlog();

  ParsedFunction<2> pf("Function", "x^2+y^2", "k=1.");

  ParameterAcceptor::initialize();

  Point<2> p(2,3);
  deallog << "F(" << p << ") = x^2+y^2+k(=1) = " << pf.value(p) << std::endl;

}
