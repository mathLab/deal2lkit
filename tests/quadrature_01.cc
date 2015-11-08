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
#include <deal2lkit/parsed_quadrature.h>
#include <deal.II/base/quadrature_lib.h>

using namespace deal2lkit;

int main ()
{
  initlog();

  ParsedQuadrature<2> pq("a", "gauss", 3);

  ParameterAcceptor::initialize();

  Quadrature<2> quadrature(pq.get_quadrature());
  QGauss<2> quad(3);

  if ( quad.size () == quadrature.size() )
    deallog << "OK! " << std::endl;

  for (unsigned int i = 0 ; i < quad.size (); ++ i)
    if ( quad.weight(i)!=quadrature.weight(i) )
      deallog << "FAIL! " << std::endl;

  deallog << "DONE! " << std::endl;
}
