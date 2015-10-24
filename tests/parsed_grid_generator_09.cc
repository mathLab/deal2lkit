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
#include <deal2lkit/parsed_grid_generator.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/base/utilities.h>


using namespace deal2lkit;

template<int dim, int spacedim>
void test(ParsedGridGenerator<dim, spacedim> &pgg)
{
  Triangulation<dim, spacedim> *tria = pgg.serial();
  GridOut go;
  go.write_msh(*tria, deallog.get_file_stream());
}


int main ()
{
  initlog();
  ParsedGridGenerator<2,2> a("Unit Hyperball", "unit_hyperball");
  ParsedGridGenerator<2,2> b("Sub Hyper Rectangle", "subhyperrectangle");
  ParsedGridGenerator<2,2> c("Hyper Shell", "hyper_shell");

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  deallog <<"Unit Hyperball"<<std::endl;
  test(a);
  deallog <<"Sub Hyper Rectangle"<<std::endl;
  test(b);
  deallog <<"Hyper Shell"<<std::endl;
  test(c);
}
