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


using namespace deal2lkit;

// Create a grid, refine it locally, write it out in ar format, read it
// back in, and check that everything is fine.

int main ()
{
  initlog();
  ParsedGridGenerator<2,2> a("Grid");

  ParameterHandler prm;
  ParameterAcceptor::declare_all_parameters(prm);

  prm.read_input_from_string(""
                             "subsection Grid\n"
                             "  set Output grid file name = grid.ar\n"
                             "end\n");

  ParameterAcceptor::parse_all_parameters(prm);

  auto t = SP(a.serial());
  t->refine_global(1);
  t->begin_active()->set_refine_flag();
  t->execute_coarsening_and_refinement();

  a.write(*t);
  std::ifstream is("grid.ar");
  std::string str((std::istreambuf_iterator<char>(is)),
                  std::istreambuf_iterator<char>());
  deallog << str;
  is.close();

  prm.read_input_from_string(""
                             "subsection Grid\n"
                             "  set Grid to generate = file\n"
                             "  set Input grid file name = grid.ar\n"
                             "  set Output grid file name = grid2.ar\n"
                             "end\n");

  ParameterAcceptor::parse_all_parameters(prm);
  auto t2 = SP(a.serial());
  a.write(*t2);
  deallog << std::endl << "========================================" << std::endl;
  is.open("grid2.ar");
  std::string str2((std::istreambuf_iterator<char>(is)),
                   std::istreambuf_iterator<char>());
  deallog << str2 << std::endl;
}
