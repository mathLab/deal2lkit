//-----------------------------------------------------------
//
//    Copyright (C) 2007 - 2015 by the deal.II authors
//    Copyright (C) 2015 - 2016 by the deal2lkit authors
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


#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal2lkit/parsed_zero_average_constraints.h>
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

  ConstraintMatrix                  cm;
  ParsedZeroAverageConstraints<dim> pnac("Parsed Zero Average Constraints",
                                         dim + 1,
                                         (dim == 2 ? "u,u,p" : "u,u,u,p"),
                                         "p");



  ParameterAcceptor::initialize();
  pnac.apply_zero_average_constraints(dof, cm);

  cm.print(deallog.get_file_stream());
}


template <int dim>
void
test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  tr.refine_global(2);

  FESystem<dim> fe(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);
  test(tr, fe);
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;
  // ParameterAcceptor::prm.log_parameters(deallog);

  test_hyper_cube<2>();
  test_hyper_cube<3>();
}
