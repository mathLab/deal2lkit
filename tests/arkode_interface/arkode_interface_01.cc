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

#include <deal2lkit/arkode_interface.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/index_set.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/numerics/error_estimator.h>

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_grid_generator.h>
#include <deal2lkit/parsed_grid_refinement.h>
#include <deal2lkit/parsed_finite_element.h>
#include <deal2lkit/error_handler.h>
#include <deal2lkit/parsed_function.h>
#include <deal2lkit/parsed_data_out.h>
#include <deal2lkit/parsed_dirichlet_bcs.h>
#include <deal2lkit/parsed_solver.h>

#include <fstream>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

using namespace dealii;
using namespace deal2lkit;

// solve heat equation

template <int dim>
class Test : public Subscriptor
{
public:

  Test();

  void run();



  shared_ptr<Vector<double> > create_new_vector() const;

  void output_step(const double t,
                   const Vector<double> &sol,
                   const unsigned int step_number);

  int explicit_rhs(const double &t,
                   const Vector<double> &y,
                   Vector<double> &expl_rhs);

  int implicit_rhs(const double &t,
                   const Vector<double> &y,
                   Vector<double> &impl_rhs);

  int assemble_mass_matrix (const double &t);

  int solve_linear_system(const double &gamma,
                          const Vector<double> &src,
                          Vector<double> &dst);

  int solve_mass_system(const Vector<double> &src,
                        Vector<double> &dst);

  int assemble_jacobian(const double &t,
                        const Vector<double> &y);

private:

  void setup_dofs();

  void make_grid();

  ConstraintMatrix constraints;

  Triangulation<dim> triangulation;
  FE_Q<dim>          fe;
  DoFHandler<dim>    dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  SparseMatrix<double> mass_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;

  ParsedFunction<dim> forcing_term;


  ARKodeInterface arkode;

};

template <int dim>
Test<dim>::Test()
  :
  fe(1),
  dof_handler(triangulation),
  forcing_term("Forcing term",1),
  arkode("ARKode parameters")
{}

template <int dim>
void Test<dim>::make_grid()
{

  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (4);
  std::cout << "   Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: "
            << triangulation.n_cells()
            << std::endl;
}

template <int dim>
void Test<dim>::setup_dofs()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit (sparsity_pattern);
  mass_matrix.reinit (sparsity_pattern);
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}

template <int dim>
int Test<dim>::assemble_mass_matrix(const double &)
{
  QGauss<dim>  quadrature_formula(2);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values   |
                           update_quadrature_points | update_JxW_values);
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_matrix = 0;
      for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_value(i, q_index) *
                                   fe_values.shape_value (j, q_index) *
                                   fe_values.JxW (q_index));
          }
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            mass_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             cell_matrix(i,j));
        }
    }

}

template <int dim>
int Test<dim>::assemble_jacobian_matrix(const double &)
{
  QGauss<dim>  quadrature_formula(2);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_gradients   |
                           update_quadrature_points | update_JxW_values);
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_matrix = 0;
      for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad (i, q_index) *
                                   fe_values.shape_grad (j, q_index) *
                                   fe_values.JxW (q_index));
          }
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));
        }
    }

}

template <int dim>
int Test<dim>::explicit_rhs(const double &t,
                            const Vector<double> &y,
                            Vector<double> &expl_rhs)
{
  expl_rhs = 0;
  forcing_term.set_time(t);
  QGauss<dim>  quadrature_formula(2);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_rhs = 0;

//      double sol;
//      Tensor<1,dim> grad_sol;

      for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
        {
//          sol = 0;
//          grad_sol=0;

//          for (unsigned int i=0; i<dofs_per_cell; ++i)
//            {
//              for (unsigned int d=0; d<dim; ++d)
//                grad_sol[d] += y[local_dof_indices[i]]*fe_values.gradient(i,q_point)[d];

//              sol += y[local_dof_indices[i]]*fe_values.shape_value(i,q_point);
//            }

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {

              cell_rhs(i) += (fe_values.shape_value (i, q_index) *
                              forcing_term.value (fe_values.quadrature_point (q_index)) *
                              fe_values.JxW (q_index));
            }
          cell->get_dof_indices (local_dof_indices);
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              expl_rhs(local_dof_indices[i]) += cell_rhs(i);
            }
        }
    }
}

template <int dim>
int Test<dim>::implicit_rhs(const double &t,
                            const Vector<double> &y,
                            Vector<double> &impl_rhs)
{
  impl_rhs = 0;
  forcing_term.set_time(t);
  QGauss<dim>  quadrature_formula(2);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_rhs = 0;

//      double sol;
      Tensor<1,dim> grad_sol;

      for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
        {
//          sol = 0;
          grad_sol=0;

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int d=0; d<dim; ++d)
                grad_sol[d] += y[local_dof_indices[i]]*fe_values.gradient(i,q_point)[d];

//              sol += y[local_dof_indices[i]]*fe_values.shape_value(i,q_point);
            }

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {

              cell_rhs(i) += (fe_values.shape_grad (i, q_index) *
                              grad_sol*
                              fe_values.JxW (q_index));
            }
          cell->get_dof_indices (local_dof_indices);
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              impl_rhs(local_dof_indices[i]) += cell_rhs(i);
            }
        }
    }
}

