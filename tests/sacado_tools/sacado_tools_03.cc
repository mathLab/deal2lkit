//-----------------------------------------------------------
//
//    Copyright (C) 2016 by the deal2lkit authors
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

// test the Val() function for SSdouble
// the .output file has been generated printing
// the variables computed directly with **Sdouble** type
// i.e., without using the Val() function.


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
#include <deal2lkit/sacado_tools.h>

#include <fstream>

#include "../tests.h"


using namespace deal2lkit;

typedef Sacado::Fad::DFad<double> Sdouble;


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
  std::vector<SSdouble> independent_local_dof_values(fe_values.dofs_per_cell);

  fe_values.reinit(dof.begin_active());
  dof.begin_active()->get_dof_indices(local_dof_indices);

  Vector<double> global_vector(dof.n_dofs());
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    global_vector[i] += i * i;

  DOFUtilities::extract_local_dofs(
    global_vector, local_dof_indices, independent_local_dof_values);

  FEValuesExtractors::Scalar scalar(dim);
  FEValuesExtractors::Vector vector(0);

  // compute divergences with sacado
  std::vector<SSdouble> div_v(quadrature.size());
  DOFUtilities::get_divergences(
    fe_values, independent_local_dof_values, vector, div_v);

  // compute gradient with sacado
  std::vector<Tensor<1, dim, SSdouble>> grad_s(quadrature.size());
  DOFUtilities::get_gradients(
    fe_values, independent_local_dof_values, scalar, grad_s);

  // compute gradient with sacado
  std::vector<Tensor<2, dim, SSdouble>> grad_v(quadrature.size());
  DOFUtilities::get_gradients(
    fe_values, independent_local_dof_values, vector, grad_v);

  // compute symmetric gradient with sacado
  std::vector<Tensor<2, dim, SSdouble>> sym_grad_v(quadrature.size());
  DOFUtilities::get_symmetric_gradients(
    fe_values, independent_local_dof_values, vector, sym_grad_v);


  {
    deallog << std::endl;
    deallog << "std::vector<Sdouble>" << std::endl;
    const std::vector<Sdouble> div_double = SacadoTools::val(div_v);
    for (unsigned int q = 0; q < quadrature.size(); ++q)
      deallog << div_double[q] << std::endl;
  }


  {
    deallog << std::endl;
    deallog << "Tensor<1," << dim << ",Sdouble>" << std::endl;
    const std::vector<Tensor<1, dim, Sdouble>> t_double =
      SacadoTools::val(grad_s);
    for (unsigned int q = 0; q < quadrature.size(); ++q)
      {
        for (unsigned int d = 0; d < dim; ++d)
          deallog << t_double[q][d] << std::endl;
      }
  }

  {
    deallog << std::endl;
    deallog << "Tensor<2," << dim << ",Sdouble>" << std::endl;
    const std::vector<Tensor<2, dim, Sdouble>> t_double =
      SacadoTools::val(grad_v);
    for (unsigned int q = 0; q < quadrature.size(); ++q)
      {
        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int dd = 0; dd < dim; ++dd)
            deallog << t_double[q][d][dd] << std::endl;
      }
  }

  {
    deallog << std::endl;
    deallog << "SymmetricTensor<2," << dim << ",Sdouble>" << std::endl;
    for (unsigned int q = 0; q < quadrature.size(); ++q)
      {
        const Tensor<2, dim, Sdouble> t_double =
          SacadoTools::val(sym_grad_v[q]);
        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int dd = 0; dd < dim; ++dd)
            deallog << t_double[d][dd] << std::endl;
      }
  }
  deallog << std::endl;
}



template <int dim>
void test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);


  FESystem<dim> fe(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);
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
