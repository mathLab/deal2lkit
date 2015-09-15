// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// check the creation of no-flux boundary conditions for a finite
// element that consists of only a single set of vector components
// (i.e. it has dim components)

#include "tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>

#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_dirichlet_bcs.h>

template<int dim>
void test (const Triangulation<dim> &tr,
           const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  deallog << "FE=" << fe.get_name()
          << std::endl;

  ConstraintMatrix cm;
//   std::set<types::boundary_id> boundary_ids;
//   for (unsigned int j=0; j<GeometryInfo<dim>::faces_per_cell; ++j)
//     boundary_ids.insert (j);
//
//   VectorTools::compute_no_normal_flux_constraints (dof, 1, boundary_ids, cm);
  ParsedDirichletBCs<dim,dim,dim+1> parsed_dirichlet("ParsedDirichletBCs",
                                                     (dim==2?"u,u,p":"u,u,u,p"),
                                                     (dim==2?"0=u.N % 1=u.N % 2=u.N % 3=u.N" :"0=u.N % 1=u.N % 2=u.N % 3=u.N % 4=u.N % 5=u.N"),
                                                     (dim==2?"0=0;0;0 % 1=1;1;0 % 2=2;2;0 % 3=3;3;0" :"0=0;0;0;0 % 1=1;1;1;0 % 2=2;2;2;0 % 3=3;3;3;0 % 4=4;4;4;0 % 5=5;5;5;0"));


  ParameterAcceptor::initialize();
  parsed_dirichlet.compute_nonzero_normal_flux_constraints(dof,cm);

  cm.print (deallog.get_file_stream ());
}


template<int dim>
void test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
    tr.begin_active()->face(i)->set_boundary_id (i);


  tr.refine_global(2);

  for (unsigned int degree=1; degree<4; ++degree)
    {
      FESystem<dim> fe (FE_Q<dim>(2), dim,
                        FE_Q<dim>(1), 1);
      test(tr, fe);
    }
}


int main()
{
  initlog();
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.threshold_double(1.e-12);
  //ParameterAcceptor::prm.log_parameters(deallog);

  test_hyper_cube<2>();
  test_hyper_cube<3>();
}
