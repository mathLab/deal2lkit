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

#include <string>
#include <fstream>
#include <streambuf>
#include <deal.II/grid/grid_out.h>


using namespace deal2lkit;

// Create a grid, refine it locally, write it out in ar format, read it
// back in, and check that everything is fine.
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
  ParsedGridGenerator<2> a("Read");

  ParameterHandler prm;
  ParameterAcceptor::declare_all_parameters(prm);
  prm.read_input_from_string(""
                             "subsection Read\n"
                             "  set Grid to generate = file \n"
                             "  set Input grid file name = "
                             SOURCE_DIR"/grids/obstacle.ucd\n"
                             "  set Manifold descriptors = 5=HyperShellBoundary\n"
                             "end\n");

  ParameterAcceptor::parse_all_parameters(prm);
  shared_ptr<Triangulation<2> > tria = SP(a.serial());
  tria->refine_global();

  GridOut grid_out;
  grid_out.write_msh (*tria, deallog.get_file_stream());
}
