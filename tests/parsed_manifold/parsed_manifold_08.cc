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
#include <deal.II/opencascade/utilities.h>
#include <deal.II/opencascade/boundary_lib.h>

#include <string>
#include <fstream>
#include <streambuf>
#include <sstream>
#include <deal.II/grid/grid_out.h>

#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax2.hxx>
#include <GC_MakeCircle.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>




using namespace deal2lkit;

// Create default manifolds for hyper shell like meshes, including
// interior parts.


using namespace OpenCASCADE;

int main ()
{
  initlog();

  std::string name = "ArclengthProjectionLineManifold";

  const unsigned int dim=2;
  const unsigned int spacedim=3;

  ParsedGridGenerator<dim,spacedim> pgg("Default");

  ParameterHandler prm;
  ParameterAcceptor::declare_all_parameters(prm);
  std::stringstream input;

  input << "subsection Default" << std::endl
        << "  set Copy boundary to manifold ids = true" << std::endl
        << "  set Colorize = true" << std::endl
        << "  set Manifold descriptors = "
        << "0=ArclengthProjectionLineManifold:"
        << SOURCE_DIR "/iges_files/edge.iges % "
        << "1=ArclengthProjectionLineManifold:"
        << SOURCE_DIR "/iges_files/edge.iges " << std::endl
        <<  "end" << std::endl;

  prm.read_input_from_string(input.str().c_str());
  ParameterAcceptor::parse_all_parameters(prm);

  shared_ptr<Triangulation<dim, spacedim> > tria = SP(pgg.serial());
  tria->refine_global(1);

  GridOut go;
  go.write_msh(*tria, deallog.get_file_stream());

  // tria->refine_global(4);
  // std::ofstream out(("/tmp/"+name+std::to_string(dim)+std::to_string(spacedim)+".msh").c_str());
  // go.write_msh(*tria, out);

  return 0;
}
