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
  std::ofstream ofile(("/tmp/mesh_"+Utilities::int_to_string(dim)
                       +Utilities::int_to_string(spacedim)+".msh").c_str());
}

int main ()
{
  initlog();
  ParsedGridGenerator<1,1> a;
  ParsedGridGenerator<1,2> b;
  ParsedGridGenerator<2,2> c;
  ParsedGridGenerator<2,3> d;
  ParsedGridGenerator<3,3> e;

  ParameterAcceptor::initialize();

  test(a);
  test(b);
  test(c);
  test(d);
  test(e);
}
