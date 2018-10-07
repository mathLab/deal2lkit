//-----------------------------------------------------------
//
//    Copyright (C) 2007 - 2015 by the deal.II authors
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


// check the creation of inhomogeneous normal-flux boundary conditions
// for a finite element that consists of only a single set of vector
// components (i.e. it has dim components)

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_dirichlet_bcs.h>
#include <deal2lkit/utilities.h>

#include <fstream>

#include "../tests.h"


using namespace deal2lkit;


template <int dim>
void
test(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  deallog << "FE=" << fe.get_name() << std::endl;

  ConstraintMatrix cm;

  // ConstantFunction<dim> constant_function(1.,dim);
  // typename FunctionMap<dim>::type function_map;
  // for (unsigned int j=0; j<GeometryInfo<dim>::faces_per_cell; ++j)
  //  function_map[j] = &constant_function;


  // std::set<types::boundary_id> boundary_ids;
  // for (unsigned int j=0; j<GeometryInfo<dim>::faces_per_cell; ++j)
  //  boundary_ids.insert (j);

  // VectorTools::compute_nonzero_normal_flux_constraints (dof, 0, boundary_ids,
  // function_map,  cm);

  ParsedDirichletBCs<dim, dim> parsed_dirichlet(
    "ParsedDirichletBCs",
    dim,
    (dim == 2 ? "u,u" : "u,u,u"),
    (dim == 2 ? "0=u.N % 1=u.N % 2=u.N % 3=u.N" :
                "0=u.N % 1=u.N % 2=u.N % 3=u.N % 4=u.N % 5=u.N"),
    (dim == 2 ? "0=1;1 % 1=1;1 % 2=1;1 % 3=1;1" :
                "0=1;1;1 % 1=1;1;1 % 2=1;1;1 % 3=1;1;1 % 4=1;1;1 % 5=1;1;1"));

  ParameterAcceptor::initialize();
  parsed_dirichlet.compute_nonzero_normal_flux_constraints(dof, cm);
  cm.print(deallog.get_file_stream());

  // Get the location of all boundary dofs
  std::vector<types::global_dof_index> face_dofs;
  const std::vector<Point<dim - 1>> &  unit_support_points =
    fe.get_unit_face_support_points();
  Quadrature<dim - 1>    quadrature(unit_support_points);
  FEFaceValues<dim, dim> fe_face_values(fe, quadrature, update_q_points);
  typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
                                                 endc = dof.end();
  for (; cell != endc; ++cell)
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
         ++face_no)
      if (cell->face(face_no)->at_boundary())
        {
          typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
          face_dofs.resize(fe.dofs_per_face);
          face->get_dof_indices(face_dofs);

          fe_face_values.reinit(cell, face_no);
          for (unsigned int i = 0; i < face_dofs.size(); ++i)
            {
              std::cout << face_dofs[i] << "\t"
                        << fe_face_values.quadrature_point(i) << "\t"
                        << fe.face_system_to_component_index(i).first
                        << std::endl;
            }
        }
}


template <int dim>
void
test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i)
    tr.begin_active()->face(i)->set_boundary_id(i);

  tr.refine_global(2);

  for (unsigned int degree = 1; degree < 4; ++degree)
    {
      FESystem<dim> fe(FE_Q<dim>(degree), dim);
      test(tr, fe);
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;
  deallog.depth_console(0);

  test_hyper_cube<2>();
  test_hyper_cube<3>();
}
