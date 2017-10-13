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


#include "../tests.h"
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
  ParsedGridGenerator<2,2> a("Rectangle");
  ParsedGridGenerator<2,2> b("Cube");

  ParameterHandler prm;
  ParameterAcceptor::declare_all_parameters(prm);

  prm.parse_input_from_string(""
                             "subsection Rectangle\n"
                             "  set Grid to generate = rectangle \n"
                             "  set Optional Point<spacedim> 1 = -1., -2. \n"
                             "  set Optional Point<spacedim> 2 =  1.,  2. \n"
                             "end\n");

  prm.parse_input_from_string(""
                             "subsection Cube\n"
                             "  set Grid to generate = rectangle \n"
                             "  set Optional Point<spacedim> 1 = -1., -1. \n"
                             "  set Optional Point<spacedim> 2 =  1.,  1. \n"
                             "end\n");

  ParameterAcceptor::parse_all_parameters(prm);

  test(a);
  test(b);
}
