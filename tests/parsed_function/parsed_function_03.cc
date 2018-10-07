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


#include <deal2lkit/parsed_function.h>

#include "../tests.h"


using namespace deal2lkit;

int
main()
{
  initlog();

  ParsedFunction<2> pf("Function", 1, "x^2+y^2");

  ParameterAcceptor::initialize();

  Point<2> p(2, 3);
  deallog << "F(" << p << ") = x^2+y^2 = " << pf.value(p) << std::endl;
}
