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

#include <fstream>

#include "../tests.h"


using namespace deal2lkit;


template <int dim>
void test(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  deallog << "FE=" << fe.get_name() << std::endl;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);


  const QGauss<dim> quadrature(2);
  FEValues<dim>     fe_values(fe, quadrature, update_values | update_gradients);

  std::vector<types::global_dof_index> local_dof_indices(
    fe_values.dofs_per_cell);
  std::vector<double> independent_local_dof_values(fe_values.dofs_per_cell);

  fe_values.reinit(dof.begin_active());
  dof.begin_active()->get_dof_indices(local_dof_indices);

  Vector<double> global_vector(dof.n_dofs());
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    global_vector[i] += i * i;

  DOFUtilities::extract_local_dofs(
    global_vector, local_dof_indices, independent_local_dof_values);

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    deallog << independent_local_dof_values[i] << std::endl;



  const unsigned int                  n_components = dim + 1;
  std::vector<std::vector<double>>    all_vars(quadrature.size(),
                                               std::vector<double>(n_components));
  std::vector<double>                 scalar_values(quadrature.size());
  std::vector<Tensor<1, dim, double>> grad_s(quadrature.size());

  std::vector<double>                 div_values(quadrature.size());
  std::vector<Tensor<1, dim, double>> vector_values(quadrature.size());
  std::vector<Tensor<2, dim, double>> grad_v(quadrature.size());
  std::vector<Tensor<2, dim, double>> sym_grad_v(quadrature.size());
  FEValuesExtractors::Scalar          scalar(dim);
  FEValuesExtractors::Vector          vector(0);

  DOFUtilities::get_values(fe_values, independent_local_dof_values, all_vars);
  DOFUtilities::get_values(
    fe_values, independent_local_dof_values, scalar, scalar_values);
  DOFUtilities::get_gradients(
    fe_values, independent_local_dof_values, scalar, grad_s);

  DOFUtilities::get_values(
    fe_values, independent_local_dof_values, vector, vector_values);
  DOFUtilities::get_divergences(
    fe_values, independent_local_dof_values, vector, div_values);
  DOFUtilities::get_gradients(
    fe_values, independent_local_dof_values, vector, grad_v);
  DOFUtilities::get_symmetric_gradients(
    fe_values, independent_local_dof_values, vector, sym_grad_v);

  deallog << "all_vars" << std::endl;
  for (unsigned int q = 0; q < quadrature.size(); ++q)
    {
      for (unsigned int i = 0; i < n_components; ++i)
        deallog << all_vars[q][i] << " ";
      deallog << std::endl;
    }

  deallog << "scalar_values" << std::endl;
  for (unsigned int q = 0; q < quadrature.size(); ++q)
    deallog << scalar_values[q] << std::endl;

  deallog << std::endl;
  deallog << "div_values" << std::endl;
  for (unsigned int q = 0; q < quadrature.size(); ++q)
    deallog << div_values[q] << std::endl;

  deallog << std::endl;
  deallog << "grad of a scalar field" << std::endl;
  for (unsigned int q = 0; q < quadrature.size(); ++q)
    {
      for (unsigned int d = 0; d < dim; ++d)
        deallog << grad_s[q][d] << " ";
      deallog << std::endl;
    }

  deallog << std::endl;
  deallog << "vector_values" << std::endl;
  for (unsigned int q = 0; q < quadrature.size(); ++q)
    {
      for (unsigned int d = 0; d < dim; ++d)
        deallog << vector_values[q][d] << " ";
      deallog << std::endl;
    }

  deallog << std::endl;
  deallog << "grad of a vector field" << std::endl;
  for (unsigned int q = 0; q < quadrature.size(); ++q)
    {
      for (unsigned int d = 0; d < dim; ++d)
        {
          for (unsigned int dd = 0; dd < dim; ++dd)
            deallog << grad_v[q][d][dd] << " ";
          deallog << std::endl;
        }
      deallog << std::endl;
    }

  deallog << std::endl;
  deallog << "sym grad of a vector field" << std::endl;
  for (unsigned int q = 0; q < quadrature.size(); ++q)
    {
      for (unsigned int d = 0; d < dim; ++d)
        {
          for (unsigned int dd = 0; dd < dim; ++dd)
            deallog << sym_grad_v[q][d][dd] << " ";
          deallog << std::endl;
        }
      deallog << std::endl;
    }
}



template <int dim>
void test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);


  FESystem<dim> fe(FE_RaviartThomas<dim>(1), 1, FE_Q<dim>(1), 1);
  test(tr, fe);
}


int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);

  deallog.attach(logfile);
  deallog.depth_console(0);

  test_hyper_cube<2>();
  test_hyper_cube<3>();
}
