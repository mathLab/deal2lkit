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

// test the DOFUtilities functions
// for Number=double

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/lac/vector.h>

#include <deal2lkit/dof_utilities.h>
#include <deal2lkit/utilities.h>

#include <fstream>

#include "../tests.h"


using namespace deal2lkit;

template <int dim>
void
test(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  deallog << "FE=" << fe.get_name() << std::endl;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);


  const QGauss<dim> quadrature(2);
  FEValues<dim>     fe_values(fe,
                          quadrature,
                          update_values | update_gradients | update_hessians);

  std::vector<types::global_dof_index> local_dof_indices(
    fe_values.dofs_per_cell);
  std::vector<double> independent_local_dof_values(fe_values.dofs_per_cell);

  fe_values.reinit(dof.begin_active());
  dof.begin_active()->get_dof_indices(local_dof_indices);

  Vector<double> global_vector(dof.n_dofs());
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    global_vector[i] += i * i;

  DOFUtilities::extract_local_dofs(global_vector,
                                   local_dof_indices,
                                   independent_local_dof_values);

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    deallog << independent_local_dof_values[i] << std::endl;

  const unsigned int                       n_components = dim + 1;
  std::vector<Tensor<1, dim>>              vector_laplacians(quadrature.size());
  std::vector<Tensor<2, dim>>              scalar_hessians(quadrature.size());
  std::vector<Tensor<3, dim>>              vector_hessians(quadrature.size());
  std::vector<Tensor<1, dim == 3 ? 3 : 1>> curls(quadrature.size());

  FEValuesExtractors::Scalar scalar(dim);
  FEValuesExtractors::Vector vector(0);

  DOFUtilities::get_laplacians(fe_values,
                               independent_local_dof_values,
                               vector,
                               vector_laplacians);
  DOFUtilities::get_hessians(fe_values,
                             independent_local_dof_values,
                             scalar,
                             scalar_hessians);
  DOFUtilities::get_hessians(fe_values,
                             independent_local_dof_values,
                             vector,
                             vector_hessians);
  DOFUtilities::get_curls(fe_values,
                          independent_local_dof_values,
                          vector,
                          curls);


  deallog << "vector_laplacians" << std::endl;
  deallog << print(vector_laplacians, "\n") << std::endl << std::endl;

  deallog << "hessians" << std::endl;
  deallog << print(scalar_hessians, "\n") << std::endl << std::endl;

  deallog << "vector_hessians" << std::endl;
  deallog << print(vector_hessians, "\n") << std::endl << std::endl;

  deallog << "curls" << std::endl;
  deallog << print(curls, "\n") << std::endl << std::endl;
}



template <int dim>
void
test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);


  FESystem<dim> fe(FE_RaviartThomas<dim>(1), 1, FE_Q<dim>(1), 1);
  test(tr, fe);


  FESystem<dim> fe2(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);

  test(tr, fe2);
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);

  deallog.attach(logfile);
  deallog.depth_console(0);

  test_hyper_cube<2>();
  test_hyper_cube<3>();
}
