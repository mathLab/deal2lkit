/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2014 by the deal.II authors
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
 * Author: Wolfgang Bangerth, Texas A&M University, 2009, 2010
 *         Timo Heister, University of Goettingen, 2009, 2010
 */

#include <deal.II/base/config.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#ifdef DEAL_II_WITH_MPI

#include <deal.II/lac/generic_linear_algebra.h>

#define USE_PETSC_LA

namespace LA
{
#ifdef USE_PETSC_LA
  using namespace dealii::LinearAlgebraPETSc;
#else
  using namespace dealii::LinearAlgebraTrilinos;
#endif
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>

#include "parsed_grid_generator.h"
#include "parsed_finite_element.h"
#include "parsed_function.h"
#include "parsed_data_out.h"
#include "utilities.h"

namespace ParallelLaplace
{
  using namespace dealii;


  template <int dim>
  class LaplaceProblem : public ParameterAcceptor
  {
  public:
    LaplaceProblem ();

    virtual void declare_parameters(ParameterHandler &prm);

    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle);
    void make_grid_fe();

    MPI_Comm                                  mpi_communicator;

    std::unique_ptr<parallel::distributed::Triangulation<dim> > triangulation;

    std::unique_ptr<FiniteElement<dim,dim> >  fe;
    std::unique_ptr<DoFHandler<dim> >         dof_handler;

    IndexSet                                  locally_owned_dofs;
    IndexSet                                  locally_relevant_dofs;

    ConstraintMatrix                          constraints;

    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector locally_relevant_solution;
    LA::MPI::Vector system_rhs;

    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;


    std::vector<unsigned int> dirichlet_boundary_ids;
    unsigned int n_cycles;
    unsigned int initial_refinement;

    // Collection of parsed_* objects
    ParsedGridGenerator<dim,dim> tria_builder;
    ParsedFiniteElement<dim,dim> fe_builder;
    ParsedFunction<dim> forcing_function;
    ParsedFunction<dim> dirichlet_function;
    ParsedDataOut<dim,dim> data_out;
  };




  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem ()
    :
    ParameterAcceptor("Global parameters"),
    mpi_communicator (MPI_COMM_WORLD),
    pcout (std::cout,
           (Utilities::MPI::this_mpi_process(mpi_communicator)
            == 0)),
    computing_timer (mpi_communicator,
                     pcout,
                     TimerOutput::summary,
                     TimerOutput::wall_times),
    tria_builder("Triangulation"),
    fe_builder("Finite element"),
    forcing_function("Rhs function"),
    dirichlet_function("Dirichlet function"),
    data_out("Data out", "vtu")
  {}


  template <int dim>
  void LaplaceProblem<dim>::declare_parameters(ParameterHandler &prm)
  {
    add_parameter(prm, &dirichlet_boundary_ids, "Dirichlet boundary ids",
                  "0", Patterns::List(Patterns::Integer(0)));

    add_parameter(prm, &n_cycles, "Number of cycles", "5",
                  Patterns::Integer(0));

    add_parameter(prm, &initial_refinement, "Initial global refinement", "3",
                  Patterns::Integer(0));

  }


  template <int dim>
  void LaplaceProblem<dim>::setup_system ()
  {
    TimerOutput::Scope t(computing_timer, "setup");

    dof_handler->distribute_dofs (*fe);

    locally_owned_dofs = dof_handler->locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (*dof_handler,
                                             locally_relevant_dofs);

    locally_relevant_solution.reinit (locally_owned_dofs,
                                      locally_relevant_dofs, mpi_communicator);
    system_rhs.reinit (locally_owned_dofs, mpi_communicator);

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (*dof_handler, constraints);
    for (auto id : dirichlet_boundary_ids)
      {
        VectorTools::interpolate_boundary_values (*dof_handler,
                                                  id,
                                                  dirichlet_function,
                                                  constraints);
      }
    constraints.close ();

    DynamicSparsityPattern csp (locally_relevant_dofs);

    DoFTools::make_sparsity_pattern (*dof_handler, csp,
                                     constraints, false);
    SparsityTools::distribute_sparsity_pattern (csp,
                                                dof_handler->n_locally_owned_dofs_per_processor(),
                                                mpi_communicator,
                                                locally_relevant_dofs);

    system_matrix.reinit (locally_owned_dofs,
                          locally_owned_dofs,
                          csp,
                          mpi_communicator);
  }




  template <int dim>
  void LaplaceProblem<dim>::assemble_system ()
  {
    TimerOutput::Scope t(computing_timer, "assembly");

    const QGauss<dim>  quadrature_formula(3);

    FEValues<dim> fe_values (*fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);

    const unsigned int   dofs_per_cell = fe->dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler->begin_active(),
    endc = dof_handler->end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          cell_rhs = 0;

          fe_values.reinit (cell);

          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            {
              const double rhs_value =
                forcing_function.value(fe_values.quadrature_point(q_point));

              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
                                         fe_values.shape_grad(j,q_point) *
                                         fe_values.JxW(q_point));

                  cell_rhs(i) += (rhs_value *
                                  fe_values.shape_value(i,q_point) *
                                  fe_values.JxW(q_point));
                }
            }

          cell->get_dof_indices (local_dof_indices);
          constraints.distribute_local_to_global (cell_matrix,
                                                  cell_rhs,
                                                  local_dof_indices,
                                                  system_matrix,
                                                  system_rhs);
        }

    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
  }




  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {
    TimerOutput::Scope t(computing_timer, "solve");
    LA::MPI::Vector
    completely_distributed_solution (locally_owned_dofs, mpi_communicator);

    SolverControl solver_control (dof_handler->n_dofs(), 1e-12);

    LA::SolverCG solver(solver_control, mpi_communicator);
    LA::MPI::PreconditionAMG preconditioner;

    LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
    data.symmetric_operator = true;
#else
    /* Trilinos defaults are good */
#endif
    preconditioner.initialize(system_matrix, data);

    solver.solve (system_matrix, completely_distributed_solution, system_rhs,
                  preconditioner);

    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;

    constraints.distribute (completely_distributed_solution);

    locally_relevant_solution = completely_distributed_solution;
  }




  template <int dim>
  void LaplaceProblem<dim>::refine_grid ()
  {
    TimerOutput::Scope t(computing_timer, "refine");

    Vector<float> estimated_error_per_cell (triangulation->n_active_cells());
    KellyErrorEstimator<dim>::estimate (*dof_handler,
                                        QGauss<dim-1>(3),
                                        typename FunctionMap<dim>::type(),
                                        locally_relevant_solution,
                                        estimated_error_per_cell);
    parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_number (*triangulation,
                                     estimated_error_per_cell,
                                     0.3, 0.03);
    triangulation->execute_coarsening_and_refinement ();
  }




  template <int dim>
  void LaplaceProblem<dim>::output_results (const unsigned int cycle)
  {
    std::stringstream suffix;
    suffix << "." << cycle;
    data_out.prepare_data_output(*dof_handler, suffix.str());
    data_out.add_data_vector (locally_relevant_solution, "u");
    data_out.write_data_and_clear();
  }

  template<int dim>
  void LaplaceProblem<dim>::make_grid_fe ()
  {
    triangulation = std::unique_ptr<parallel::distributed::Triangulation<dim> >
                    (tria_builder.distributed(mpi_communicator));
    dof_handler = std::unique_ptr<DoFHandler<dim> >
                  (new DoFHandler<dim>(*triangulation));

    triangulation->refine_global(initial_refinement);

    std::cout << "Number of active cells: "
              << triangulation->n_active_cells()
              << std::endl;

    fe=std::unique_ptr<FiniteElement<dim> >
       (fe_builder());
  }


  template <int dim>
  void LaplaceProblem<dim>::run ()
  {
    make_grid_fe();
    for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
      {
        pcout << "Cycle " << cycle << ':' << std::endl;

        if (cycle != 0)
          refine_grid ();

        setup_system ();

        pcout << "   Number of active cells:       "
              << triangulation->n_global_active_cells()
              << std::endl
              << "   Number of degrees of freedom: "
              << dof_handler->n_dofs()
              << std::endl;

        assemble_system ();
        solve ();

        if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
          {
            TimerOutput::Scope t(computing_timer, "output");
            output_results (cycle);
          }

        computing_timer.print_summary ();
        computing_timer.reset ();

        pcout << std::endl;
      }
  }
}




int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace ParallelLaplace;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      deallog.depth_console (0);

      {
        LaplaceProblem<2> laplace_problem_2d;
        ParameterAcceptor::initialize("parameters.prm", "used_parameters.prm");
        laplace_problem_2d.run ();
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}

#else
int main()
{
  return 0;
}
#endif
