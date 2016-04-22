//-----------------------------------------------------------
//
//    Copyright (C) 2016 by the deal2lkit authors
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

// test the DOFUtilities functions
// for Number=double

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_parser.h>
#include <deal2lkit/dof_utilities.h>


using namespace deal2lkit;


int main()
{
  std::ofstream logfile ("output");

  deallog.attach(logfile);
  deallog.depth_console (0);

  // set up problem:
  std::string variables = "t";
  std::string _step_size = "(t<1?1e-3:(t<2?2e-3:7))";
  std::map<std::string,double> constants;
  // FunctionParser with 2 variables and 1 component:
  FunctionParser<1> fp(1);
  fp.initialize(variables,
                _step_size,
                constants);
  // Point at which we want to evaluate the function
  Point<1> point1(0.0);
  // evaluate the expression at 'point':
  double result = fp.value(point1);
  deallog << result << std::endl;

  Point<1> point2(1.0);
  result = fp.value(point2);
  deallog << result << std::endl;

  Point<1> point3(56);
  result = fp.value(point3);
  deallog << result << std::endl;
}
