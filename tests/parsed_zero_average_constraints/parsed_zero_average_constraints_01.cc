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

// test basic functionalities


#include <deal2lkit/parsed_zero_average_constraints.h>

#include "../tests.h"


using namespace deal2lkit;


int
main()
{
  initlog();
  ParsedZeroAverageConstraints<2, 2> pnac("ciao", 1, "u", "u", "");

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);
}
