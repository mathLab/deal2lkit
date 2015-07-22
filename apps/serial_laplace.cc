/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2014 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 */



#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>


#include <deal.II/base/mpi.h>

#include <fstream>
#include <iostream>


#include "parsed_grid_generator.h"
#include "parsed_finite_element.h"
#include "parsed_dirichlet_bcs.h"
#include "parsed_function.h"
#include "parsed_data_out.h"
#include "error_handler.h"
#include "utilities.h"

using namespace dealii;



template<int dim>
class SerialLaplace : public ParameterAcceptor
{
public:

  SerialLaplace ();

  virtual void declare_parameters(ParameterHandler &prm);

  void run ();


private:
  void make_grid_fe ();
  void setup_system ();
  void assemble_system ();
  void solve ();
  void output_results ();
  void process_solution (const unsigned int cycle);
  void refine_grid ();

  shared_ptr<Triangulation<dim> >    triangulation;
  shared_ptr<FiniteElement<dim,dim> >             fe;
  shared_ptr<DoFHandler<dim> >       dof_handler;

  ConstraintMatrix     constraints;
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;

  ParsedGridGenerator<dim,dim> tria_builder;
  ParsedFiniteElement<dim,dim> fe_builder;

  ParsedFunction<dim> forcing_function;
  ParsedFunction<dim> permeability;
  ParsedFunction<dim> exact_solution;
  ParsedDirichletBCs<dim,dim,1> dirichlet_bcs;

  ParsedDataOut<dim,dim> data_out;

  unsigned int n_cycles;
  unsigned int initial_refinement;
  std::string   refinement_mode;
  ErrorHandler<2> eh;
};


template <int dim>
SerialLaplace<dim>::SerialLaplace ()
  :
  ParameterAcceptor("Global parameters"),
  tria_builder("Triangulation"),
  fe_builder("Finite element"),
  forcing_function("Rhs function","2 * pi * pi * sin(pi*x) * sin(pi*y)"),
  permeability("Permeability","1."),
  exact_solution("Exact solution","sin(pi*x) * sin(pi*y)"),
  dirichlet_bcs("Dirichlet BCs", "u", "0=u", "0=0"),
  data_out("Data out", "vtu")
{}


template <int dim>
void SerialLaplace<dim>::declare_parameters(ParameterHandler &prm)
{

  add_parameter(prm, &initial_refinement, "Initial global refinement", "4",
                Patterns::Integer(0));

  add_parameter(prm, &n_cycles, "Total number of cycles", "5",
                Patterns::Integer(0));

  add_parameter(prm, &refinement_mode, "Refinement strategy", "adaptive_refinement",
                Patterns::Selection("adaptive_refinement|global_refinement"));

}

template <int dim>
void SerialLaplace<dim>::make_grid_fe ()
{


  triangulation = SP(tria_builder.serial());

  triangulation->refine_global(initial_refinement);

  dof_handler = SP(new DoFHandler<dim>(*triangulation));

  std::cout << "Number of active cells: "
            << triangulation->n_active_cells()
            << std::endl;

  fe=SP(fe_builder());

}



template <int dim>
void SerialLaplace<dim>::setup_system ()
{
  dof_handler->distribute_dofs (*fe);
  std::cout << "Number of degrees of freedom: "
            << dof_handler->n_dofs()
            << std::endl;

  constraints.clear ();

  DoFTools::make_hanging_node_constraints (*dof_handler, constraints);

  dirichlet_bcs.interpolate_boundary_values(*dof_handler,constraints);

  constraints.close ();

  DynamicSparsityPattern d_sparsity(dof_handler->n_dofs());

  DoFTools::make_sparsity_pattern (*dof_handler, d_sparsity,
                                   constraints);

  sparsity_pattern.copy_from(d_sparsity);

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler->n_dofs());

  system_rhs.reinit (dof_handler->n_dofs());
}


template <int dim>
void SerialLaplace<dim>::assemble_system ()
{
  QGauss<dim>  quadrature_formula(4);

  MatrixCreator::create_laplace_matrix (StaticMappingQ1<dim, dim>::mapping,
                                        *dof_handler, quadrature_formula, system_matrix,
                                        forcing_function, system_rhs, &permeability, constraints);

}


template <int dim>
void SerialLaplace<dim>::solve ()
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              solver (solver_control);
  solver.solve (system_matrix, solution, system_rhs,
                PreconditionIdentity());
  constraints.distribute (solution);
}


template <int dim>
void SerialLaplace<dim>::output_results ()
{
  data_out.prepare_data_output(*dof_handler);
  data_out.add_data_vector (solution, "u");
  data_out.write_data_and_clear();
}

template <int dim>
void SerialLaplace<dim>::process_solution (const unsigned int cycle)
{
  eh.error_from_exact(*dof_handler, solution, exact_solution);
}

template <int dim>
void SerialLaplace<dim>::refine_grid ()
{
  if (refinement_mode == "global_refinement")
    {
      triangulation->refine_global (1);
    }
  else if (refinement_mode == "adaptive_refinement")
    {
      Vector<float> estimated_error_per_cell (triangulation->n_active_cells());

      KellyErrorEstimator<dim>::estimate (*dof_handler,
                                          QGauss<dim-1>(3),
                                          typename FunctionMap<dim>::type(),
                                          solution,
                                          estimated_error_per_cell);

      GridRefinement::refine_and_coarsen_fixed_number (*triangulation,
                                                       estimated_error_per_cell,
                                                       0.3, 0.03);

      triangulation->execute_coarsening_and_refinement ();
    }
  else
    {
      Assert (false, ExcNotImplemented());
    }
}

template <int dim>
void SerialLaplace<dim>::run ()
{
  make_grid_fe ();
  for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
    {
      std::cout<<"cycle : "<<cycle+1<<std::endl;
      setup_system ();
      assemble_system ();
      solve ();
      process_solution (cycle);
      if (cycle < n_cycles-1)
        refine_grid ();
      else
        output_results ();
    }
  eh.output_table(std::cout);

}



int main (int argc, char *argv[])
{

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  SerialLaplace<2> laplace_problem_2;
  ParameterAcceptor::initialize("parameters_ser.prm", "used_parameters_ser.prm");
  laplace_problem_2.run ();

  return 0;
}
