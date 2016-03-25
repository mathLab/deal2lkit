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

// test the inner functions

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/base/exceptions.h>

#include <fstream>

#include <deal2lkit/dof_utilities.h>
#include <deal2lkit/sacado_tools.h>

using namespace deal2lkit;

typedef Sacado::Fad::DFad<double> Sdouble;
typedef Sacado::Fad::DFad<Sdouble> SSdouble;


template<int dim>
void test (const Triangulation<dim> &tr,
           const FiniteElement<dim> &fe)
{
  deallog << "FE=" << fe.get_name()
          << std::endl;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);


  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (fe, quadrature,
                           update_values | update_gradients);

  std::vector<types::global_dof_index>    local_dof_indices (fe_values.dofs_per_cell);
  std::vector<double> ildv_double (fe_values.dofs_per_cell);
  std::vector<Sdouble> ildv_sdouble (fe_values.dofs_per_cell);
  std::vector<SSdouble> ildv_ssdouble (fe_values.dofs_per_cell);

  fe_values.reinit (dof.begin_active());
  dof.begin_active()->get_dof_indices (local_dof_indices);

  Vector<double> global_vector(dof.n_dofs());
  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    global_vector[i] += i*i;

  DOFUtilities::extract_local_dofs(global_vector, local_dof_indices, ildv_double);
  DOFUtilities::extract_local_dofs(global_vector, local_dof_indices, ildv_sdouble);
  DOFUtilities::extract_local_dofs(global_vector, local_dof_indices, ildv_ssdouble);


  FEValuesExtractors::Scalar scalar (dim);
  FEValuesExtractors::Vector vector (0);

// double
  std::vector <Tensor <1, dim, double> > grad_s_double(quadrature.size());
  std::vector <Tensor <1, dim, double> > vector_values_double(quadrature.size());
  std::vector <Tensor <2, dim, double> > grad_v_double(quadrature.size());
  std::vector <Tensor <2, dim, double> > sym_grad_v_double(quadrature.size());

  DOFUtilities::get_gradients(fe_values, ildv_double, scalar, grad_s_double);
  DOFUtilities::get_values(fe_values, ildv_double, vector, vector_values_double);
  DOFUtilities::get_gradients(fe_values, ildv_double, vector, grad_v_double);
  DOFUtilities::get_symmetric_gradients(fe_values, ildv_double, vector, sym_grad_v_double);

// Sdouble
  std::vector <Tensor <1, dim, Sdouble> > grad_s_sdouble(quadrature.size());
  std::vector <Tensor <1, dim, Sdouble> > vector_values_sdouble(quadrature.size());
  std::vector <Tensor <2, dim, Sdouble> > grad_v_sdouble(quadrature.size());
  std::vector <Tensor <2, dim, Sdouble> > sym_grad_v_sdouble(quadrature.size());

  DOFUtilities::get_gradients(fe_values, ildv_sdouble, scalar, grad_s_sdouble);
  DOFUtilities::get_values(fe_values, ildv_sdouble, vector, vector_values_sdouble);
  DOFUtilities::get_gradients(fe_values, ildv_sdouble, vector, grad_v_sdouble);
  DOFUtilities::get_symmetric_gradients(fe_values, ildv_sdouble, vector, sym_grad_v_sdouble);

// SSdouble
  std::vector <Tensor <1, dim, SSdouble> > grad_s_ssdouble(quadrature.size());
  std::vector <Tensor <1, dim, SSdouble> > vector_values_ssdouble(quadrature.size());
  std::vector <Tensor <2, dim, SSdouble> > grad_v_ssdouble(quadrature.size());
  std::vector <Tensor <2, dim, SSdouble> > sym_grad_v_ssdouble(quadrature.size());

  DOFUtilities::get_gradients(fe_values, ildv_ssdouble, scalar, grad_s_ssdouble);
  DOFUtilities::get_values(fe_values, ildv_ssdouble, vector, vector_values_ssdouble);
  DOFUtilities::get_gradients(fe_values, ildv_ssdouble, vector, grad_v_ssdouble);
  DOFUtilities::get_symmetric_gradients(fe_values, ildv_ssdouble, vector, sym_grad_v_ssdouble);

  std::vector<double> res_double (quadrature.size(),0);
  std::vector<Sdouble> res_sdouble (quadrature.size(),0);
  std::vector<SSdouble> res_ssdouble (quadrature.size(),0);

  for (unsigned int q=0; q<quadrature.size(); ++q)
    {
      res_double[q] = grad_s_double[q]*grad_s_double[q];
      deallog << res_double[q] <<std::endl;
      res_sdouble[q] = SacadoTools::scalar_product(grad_s_double[q],grad_s_sdouble[q]);
      deallog << res_sdouble[q] <<std::endl;
      res_ssdouble[q] = SacadoTools::scalar_product(grad_s_ssdouble[q],grad_s_double[q]);
      deallog << res_ssdouble[q] <<std::endl;
      deallog << res_double[q] - res_sdouble[q].val() <<std::endl;
      deallog << res_double[q] - res_ssdouble[q].val().val() <<std::endl;
    }

  for (unsigned int q=0; q<quadrature.size(); ++q)
    {
      res_double[q] = vector_values_double[q]*vector_values_double[q];
      deallog << res_double[q] <<std::endl;
      res_sdouble[q] = SacadoTools::scalar_product(vector_values_double[q],vector_values_sdouble[q]);
      deallog << res_sdouble[q] <<std::endl;
      res_ssdouble[q] = SacadoTools::scalar_product(vector_values_ssdouble[q],vector_values_double[q]);
      deallog << res_ssdouble[q] <<std::endl;
      deallog << res_double[q] - res_sdouble[q].val() <<std::endl;
      deallog << res_double[q] - res_ssdouble[q].val().val() <<std::endl;
    }

  for (unsigned int q=0; q<quadrature.size(); ++q)
    {
      res_double[q] = scalar_product(grad_v_double[q],grad_v_double[q]);
      deallog << res_double[q] <<std::endl;
      res_sdouble[q] = SacadoTools::scalar_product(grad_v_sdouble[q],grad_v_double[q]);
      deallog << res_sdouble[q] <<std::endl;
      res_ssdouble[q] = SacadoTools::scalar_product(grad_v_double[q],grad_v_ssdouble[q]);
      deallog << res_ssdouble[q] <<std::endl;
      deallog << res_double[q] - res_sdouble[q].val() <<std::endl;
      deallog << res_double[q] - res_ssdouble[q].val().val() <<std::endl;
    }

  for (unsigned int q=0; q<quadrature.size(); ++q)
    {
      res_double[q] = scalar_product(sym_grad_v_double[q],sym_grad_v_double[q]);
      deallog << res_double[q] <<std::endl;
      res_sdouble[q] = SacadoTools::scalar_product(sym_grad_v_sdouble[q],sym_grad_v_double[q]);
      deallog << res_sdouble[q] <<std::endl;
      res_ssdouble[q] = SacadoTools::scalar_product(sym_grad_v_double[q],sym_grad_v_ssdouble[q]);
      deallog << res_ssdouble[q] <<std::endl;
      deallog << res_double[q] - res_sdouble[q].val() <<std::endl;
      deallog << res_double[q] - res_ssdouble[q].val().val() <<std::endl;
    }

  for (unsigned int q=0; q<quadrature.size(); ++q)
    {
      res_double[q] = SacadoTools::scalar_product(vector_values_double[q],vector_values_double[q]);
      res_double[q] = SacadoTools::scalar_product(grad_v_double[q],grad_v_double[q]);
    }


}



template<int dim>
void test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);


  FESystem<dim> fe (FE_Q<dim>(2), dim,
                    FE_Q<dim>(1), 1);
  test(tr, fe);
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test_hyper_cube<2>();
  test_hyper_cube<3>();
}
