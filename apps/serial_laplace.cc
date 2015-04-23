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

#include <deal.II/fe/fe_values.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
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
#include "parsed_function.h"
#include "parsed_data_out.h"
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
  ParsedFunction<dim> dirichlet_function;

  ParsedDataOut<dim,dim> data_out;

  std::vector<unsigned int> dirichlet_boundary_ids;
  unsigned int n_cycles;
  unsigned int initial_refinement;
};


template <int dim>
SerialLaplace<dim>::SerialLaplace ()
  :
  ParameterAcceptor("Global parameters"),
  tria_builder("Triangulation"),
  fe_builder("Finite element"),
  forcing_function("Rhs function"),
  permeability("Permeability","1."),
  dirichlet_function("Dirichlet function"),
  data_out("Data out", "vtu")
{}


template <int dim>
void SerialLaplace<dim>::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &dirichlet_boundary_ids, "Dirichlet boundary ids",
                "0", Patterns::List(Patterns::Integer(0)));

  add_parameter(prm, &initial_refinement, "Initial global refinement", "4",
                Patterns::Integer(0));

}

template <int dim>
void SerialLaplace<dim>::make_grid_fe ()
{

  ParameterAcceptor::initialize("parameters_ser.prm", "used_parameters_ser.prm");

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

  for (auto id : dirichlet_boundary_ids)
    {
      VectorTools::interpolate_boundary_values (*dof_handler,
                                                id,
                                                dirichlet_function,
                                                constraints);
    }
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
}


template <int dim>
void SerialLaplace<dim>::output_results ()
{
  data_out.prepare_data_output(*dof_handler);
  data_out.add_data_vector (solution, "u");
  data_out.write_data_and_clear();
}


template <int dim>
void SerialLaplace<dim>::run ()
{
  make_grid_fe ();
  setup_system ();
  assemble_system ();
  solve ();
  output_results ();
}



int main (int argc, char *argv[])
{

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  SerialLaplace<2> laplace_problem_2;
  laplace_problem_2.run ();

  return 0;
}
