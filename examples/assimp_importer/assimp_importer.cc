/* ---------------------------------------------------------------------
 * Copyright (C) 2015 by the deal2lkit authors
 *
 * This file is part of the deal2lkit library.
 * ---------------------------------------------------------------------
 */

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>

#include <deal2lkit/config.h>
#include <deal2lkit/assimp_interface.h>

#include <fstream>

using namespace dealii;

// A minimal program that takes as input a file in a format supported
// by Assimp, and generates a grid supported by deal.II

int main (int argc, char ** argv)
{
  AssertThrow(argc == 3, ExcMessage("Must specify exactly 2 arguments!"));
  Triangulation<2,3> tria;
  if(AssimpInterface::generate_triangulation(argv[1], tria, -1, true, 1e-3)) {
    GridOut grid_out;
    std::string fname(argv[2]);
    std::ofstream out (argv[2]);
    
    AssertThrow(out, ExcIO());
    if(fname.find("msh") !=  std::string::npos) {
      grid_out.write_msh (tria, out);
    } else if(fname.find("vtk") !=  std::string::npos) {
      grid_out.write_vtk (tria, out);
    } else if(fname.find("inp") != std::string::npos) {
      grid_out.write_ucd (tria, out);
    } else {
      AssertThrow(false, ExcMessage("Unrecognized output format. vtk, msh or inp."));
    }
  }
}
