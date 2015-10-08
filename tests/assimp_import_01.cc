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
#include <deal.II/grid/tria.h>

#include <deal2lkit/assimp_interface.h>

int main ()
{
  Triangulation<2,3> tria;
  AssimpInterface::generate_triangulation(SOURCE_DIR "/grids/torus.obj", tria, -1, true, 1e-3);
  return 0;
}
