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
#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_finite_element.h>


int main ()
{
  initlog();
  ParsedFiniteElement<2> fe_builder("FE", "FESystem[FE_Q(2)^d-FE_DGP(1)]", "u,u,p", 1);

  ParameterAcceptor::initialize("test_in.prm", "test_out.prm");
  // Should fail because fe_name has 3 components
  return 0;
}
