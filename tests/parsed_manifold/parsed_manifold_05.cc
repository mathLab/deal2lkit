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

#include <string>
#include <fstream>
#include <streambuf>
#include <sstream>
#include <deal.II/grid/grid_out.h>


using namespace deal2lkit;

// Create a grid, refine it locally, write it out in ar format, read it
// back in, and check that everything is fine.
template<int dim, int spacedim>
void test(const std::string &name)
{
  deallog << "Testing " << name
          <<"<"<<dim<<","<<spacedim<<">"<< std::endl;

  ParsedGridGenerator<dim, spacedim> pgg("Default");

  ParameterHandler prm;
  ParameterAcceptor::declare_all_parameters(prm);
  std::stringstream input;
  unsigned int ncells = name=="half_hyper_shell" ? 4 : 6;
  if (name == "cylinder")
    ncells = 0;
  if (name == "cylinder_shell")
    ncells = 2;
  input <<  "subsection Default" << std::endl
        <<  "  set Grid to generate = " << name << std::endl
        <<  "  set Colorize = false" << std::endl
        <<  "  set Create default manifolds = true" << std::endl
        <<  "  set Optional int 1 = " << ncells << std::endl
        <<  "  set Optional double 3 = 1.0" << std::endl
        <<  "  set Optional Point<spacedim> 1= -1,0,0" << std::endl
        <<  "  set Optional Point<spacedim> 2=  1,0,0" << std::endl
        <<  "end" << std::endl;
  prm.read_input_from_string(input.str().c_str());
  ParameterAcceptor::parse_all_parameters(prm);

  shared_ptr<Triangulation<dim, spacedim> > tria = SP(pgg.serial());
  tria->refine_global(1);

  GridOut go;
  go.write_msh(*tria, deallog.get_file_stream());
}

int main ()
{
  initlog();

  test<3,3>("truncated_cone");

}
