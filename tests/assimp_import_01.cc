//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

#include "tests.h"
#include <deal.II/grid/tria.h>

#include <deal2lkit/assimp_interface.h>

int main ()
{
  Triangulation<2,3> tria;
  AssimpInterface::generate_triangulation(SOURCE_DIR "/grids/torus.obj", tria, -1, true, 1e-3);
  return 0;
}
