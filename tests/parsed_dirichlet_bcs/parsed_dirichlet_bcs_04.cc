//-----------------------------------------------------------
//
//    Copyright (C) 2006 - 2014 by the deal2lkit authors
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


// VectorTools::interpolate_boundary_values still had bugs in 1d after
// switching to a scheme where we can assign boundary indicators also in 1d

#include <deal.II/base/function_lib.h>
#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_dirichlet_bcs.h>
#include <deal2lkit/utilities.h>

#include <fstream>
#include <iomanip>
#include <vector>

#include "../tests.h"


using namespace deal2lkit;


template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;

  Triangulation<dim> triangulation;
  FE_Q<dim>          fe(1);
  DoFHandler<dim>    dof_handler(triangulation);

  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.begin_active()->face(0)->set_boundary_id(10);
  triangulation.begin_active()->face(1)->set_boundary_id(20);
  triangulation.refine_global(1);

  dof_handler.distribute_dofs(fe);
  QGauss<dim - 1> quadrature(2);
  MappingQ1<dim>  mapping;

  typename FunctionMap<dim>::type boundary_map;
  Functions::SquareFunction<dim>  f;
  boundary_map[10] = &f;
  boundary_map[20] = &f;


  ParsedDirichletBCs<dim, dim> parsed_dirichlet(
    "Parsed Dirichlet BCs",
    1,
    "",
    "10=0 % 20=0",
    (dim == 1 ? "10=x^2 % 20=x^2" :
                (dim == 2 ? "10=x^2+y^2 % 20=x^2+y^2" :
                            "10=x^2+y^2+z^2 % 20=x^2+y^2+z^2")));

  dealii::ParameterAcceptor::initialize();
  std::map<types::global_dof_index, double> boundary_values;
  parsed_dirichlet.project_boundary_values(mapping,
                                           dof_handler,
                                           quadrature,
                                           boundary_values);
  //  VectorTools::project_boundary_values(mapping, dof_handler, boundary_map,
  //  quadrature,
  //                                       boundary_values);
  deallog << boundary_values.size() << std::endl;
  for (std::map<types::global_dof_index, double>::const_iterator p =
         boundary_values.begin();
       p != boundary_values.end();
       ++p)
    deallog << p->first << ' ' << p->second << std::endl;
}


int
main()
{
  initlog();
  deallog.depth_console(0);
  dealii::ParameterAcceptor::prm.log_parameters(deallog);

  test<1>();
  test<2>();
  test<3>();
}
