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

#include "heat_ida.h"

#ifdef D2K_WITH_SUNDIALS

#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <nvector/nvector_parallel.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/packaged_operation.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/distributed/solution_transfer.h>
//#include <deal.II/base/index_set.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_math.h>

template <int dim>
Heat<dim>::Heat (const MPI_Comm &communicator)
  :
  SundialsInterface<VEC> (communicator),
  comm(communicator),
  pcout (std::cout,
         (Utilities::MPI::this_mpi_process(comm)
          == 0)),
  timer_outfile("timer.txt"),
  tcout (timer_outfile,
         (Utilities::MPI::this_mpi_process(comm)
          == 0)),
  computing_timer (comm,
                   tcout,
                   TimerOutput::summary,
                   TimerOutput::wall_times),

  eh("Error Tables", "u",
     "L2,H1"),

  pgg("Domain"),
  pgr("Refinement"),
  fe_builder("Finite Element"),

  exact_solution("Exact solution"),
  forcing_term("Forcing term"),
  initial_solution("Initial solution"),
  initial_solution_dot("Initial solution_dot"),
  dirichlet_bcs("Dirichlet BCs", "u", "0=u"),

  data_out("Output Parameters", "vtu"),
  Ainv( "Solver", "cg",
      /* iter= */ 1000,
      /* reduction= */1e-8,
        linear_operator<VEC>(jacobian_matrix) ),
  dae(*this)
{}

template <int dim>
void Heat<dim>::declare_parameters (ParameterHandler &prm)
{
  add_parameter(  prm,
                  &initial_global_refinement,
                  "Initial global refinement",
                  "1",
                  Patterns::Integer (0));

  add_parameter(  prm,
                  &n_cycles,
                  "Number of cycles",
                  "3",
                  Patterns::Integer (0));

  add_parameter(  prm,
                  &max_time_iterations,
                  "Maximum number of time steps",
                  "10000",
                  Patterns::Integer (0));

  add_parameter(  prm,
                  &timer_file_name,
                  "Timer output file",
                  "timer.txt",
                  Patterns::FileName());


  add_parameter(  prm,
                  &adaptive_refinement,
                  "Adaptive refinement",
                  "true",
                  Patterns::Bool());


  add_parameter(  prm,
                  &use_direct_solver,
                  "Use direct solver if available",
                  "true",
                  Patterns::Bool());

  add_parameter(  prm,
                  &use_space_adaptivity,
                  "Refine mesh during transient",
                  "true",
                  Patterns::Bool());

  add_parameter(  prm,
                  &kelly_threshold,
                  "Threshold for restart solver",
                  "1e-2",
                  Patterns::Double(0.0));

  add_parameter(  prm,
                  &max_cells,
                  "Maximum number of cells",
                  "1000",
                  Patterns::Integer(),
                  "If negative, there is no upper bound");

  add_parameter(  prm,
                  &top_fraction,
                  "Top fraction",
                  "0.3",
                  Patterns::Double(0.0));

  add_parameter(  prm,
                  &bottom_fraction,
                  "Bottom fraction",
                  "0.1",
                  Patterns::Double(0.0));

  add_parameter(  prm,
                  &diffusivity,
                  "Diffusivity",
                  "1.0",
                  Patterns::Double(0.0));
}

template <int dim>
void Heat<dim>::make_grid_fe()
{
  triangulation = SP(pgg.distributed(comm));
  dof_handler = SP(new DoFHandler<dim>(*triangulation));
  fe=SP(fe_builder());
  triangulation->refine_global (initial_global_refinement);
}


template <int dim>
void Heat<dim>::setup_dofs (const bool &first_run)
{
  computing_timer.enter_section("Setup dof systems");

  dof_handler->distribute_dofs (*fe);

  mapping = SP(new MappingQ<dim>(1));

  const unsigned int n_dofs = dof_handler->n_dofs();

  std::locale s = pcout.get_stream().getloc();
  pcout.get_stream().imbue(std::locale(""));
  pcout << "Number of active cells: "
        << triangulation->n_global_active_cells()
        << " (on "
        << triangulation->n_levels()
        << " levels)"
        << std::endl
        << "Number of degrees of freedom: "
        << n_dofs
        << std::endl;
  pcout.get_stream().imbue(s);


  partitioning = dof_handler->locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(*dof_handler,
                                          relevant_partitioning);

  constraints.clear ();
  constraints.reinit (relevant_partitioning);

  DoFTools::make_hanging_node_constraints (*dof_handler,
                                           constraints);

  dirichlet_bcs.interpolate_boundary_values(*dof_handler, constraints);
  constraints.close ();

  jacobian_matrix.clear();
  jacobian_matrix_sp.reinit(partitioning,partitioning,relevant_partitioning,comm);


  DoFTools::make_sparsity_pattern (*dof_handler,
                                   jacobian_matrix_sp,
                                   constraints,
                                   false,
                                   Utilities::MPI::this_mpi_process(comm));

  jacobian_matrix_sp.compress();

  jacobian_matrix.reinit(jacobian_matrix_sp);

  solution.reinit(partitioning, comm);
  solution_dot.reinit(partitioning, comm);

  distributed_solution.reinit(partitioning,relevant_partitioning,comm);
  distributed_solution_dot.reinit(partitioning,relevant_partitioning,comm);

  if (first_run)
    {
      VectorTools::interpolate(*dof_handler, initial_solution, solution);
      VectorTools::interpolate(*dof_handler, initial_solution_dot, solution_dot);
    }

  computing_timer.exit_section();
}


template <int dim>
void Heat<dim>::assemble_jacobian_matrix(const double t,
                                         const VEC &solution,
                                         const VEC &solution_dot,
                                         const double alpha)
{

  computing_timer.enter_section ("   Assemble jacobian matrix");
  jacobian_matrix = 0;
  dirichlet_bcs.set_time(t);
  exact_solution.set_time(t);
  constraints.clear();
  DoFTools::make_hanging_node_constraints (*dof_handler,
                                           constraints);

  dirichlet_bcs.interpolate_boundary_values(*dof_handler, constraints);

  constraints.close ();

  VEC tmp(solution);
  constraints.distribute(tmp);
  distributed_solution = tmp;
  distributed_solution_dot = solution_dot;

  const QGauss<dim>  quadrature_formula(fe->degree+1);

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

        fe_values.reinit (cell);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          {

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  cell_matrix(i,j) += (
                                        alpha*fe_values.shape_value(i,q_point) *
                                        fe_values.shape_value(j,q_point)

                                        +
                                        diffusivity *
                                        fe_values.shape_grad(i,q_point) *
                                        fe_values.shape_grad(j,q_point)

                                      )*fe_values.JxW(q_point);

              }
          }

        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (cell_matrix,
                                                local_dof_indices,
                                                jacobian_matrix);
      }

  jacobian_matrix.compress (VectorOperation::add);

  auto id = solution.locally_owned_elements();
  for (unsigned int i=0; i<id.n_elements(); ++i)
    {
      auto j = id.nth_index_in_set(i);
      if (constraints.is_constrained(j))
        jacobian_matrix.set(j, j, 1.0);
    }
  //compress(jacobian_matrix,VectorOperation::insert);
  jacobian_matrix.compress(VectorOperation::insert);

  computing_timer.exit_section();
}

template <int dim>
int Heat<dim>::residual (const double t,
                         const VEC &solution,
                         const VEC &solution_dot,
                         VEC &dst)
{
  computing_timer.enter_section ("Residual");
  dirichlet_bcs.set_time(t);
  forcing_term.set_time(t);
  exact_solution.set_time(t);
  constraints.clear();
  DoFTools::make_hanging_node_constraints (*dof_handler,
                                           constraints);

  dirichlet_bcs.interpolate_boundary_values(*dof_handler, constraints);

  constraints.close ();

  VEC tmp(solution);
  constraints.distribute(tmp);

  distributed_solution = tmp;
  distributed_solution_dot = solution_dot;

  dst = 0;

  const QGauss<dim>  quadrature_formula(fe->degree+1);

  FEValues<dim> fe_values (*fe, quadrature_formula,
                           update_values |  update_gradients |
                           update_quadrature_points |
                           update_JxW_values);

  const unsigned int   dofs_per_cell = fe->dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  const FEValuesExtractors::Scalar u (0);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler->begin_active(),
  endc = dof_handler->end();
  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned())
      {
        cell_rhs = 0;

        fe_values.reinit (cell);
        cell->get_dof_indices (local_dof_indices);

        std::vector<Point<dim> > quad_points(n_q_points);

        quad_points = fe_values.get_quadrature_points();

        double sol_dot;
        Tensor<1,dim> grad_sol;


        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          {
            grad_sol = 0.0;
            sol_dot = 0.0;
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int d=0; d<dim; ++d)
                  grad_sol[d] += distributed_solution[local_dof_indices[i]]*fe_values[u].gradient(i,q_point)[d];

                sol_dot += distributed_solution_dot[local_dof_indices[i]]*fe_values.shape_value(i,q_point);
              }

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {

                cell_rhs(i) += (
                                 sol_dot*
                                 fe_values.shape_value(i,q_point)

                                 +
                                 diffusivity *
                                 grad_sol *
                                 fe_values.shape_grad(i,q_point)

                                 -
                                 forcing_term.value(quad_points[q_point]) *
                                 fe_values.shape_value(i,q_point)

                               )*fe_values.JxW(q_point);
              }
          }

        constraints.distribute_local_to_global (cell_rhs,
                                                local_dof_indices,
                                                dst);
      }

  dst.compress (VectorOperation::add);

  auto id = solution.locally_owned_elements();
  for (unsigned int i=0; i<id.n_elements(); ++i)
    {
      auto j = id.nth_index_in_set(i);
      if (constraints.is_constrained(j))
        dst[j] = solution(j)-distributed_solution(j);
    }

  dst.compress(VectorOperation::insert);
  computing_timer.exit_section();
  return 0;
}

template <int dim>
shared_ptr<VEC> Heat<dim>::create_new_vector () const
{
  shared_ptr<VEC> ret = SP(new VEC(solution));
  return ret;
}

template <int dim>
unsigned int Heat<dim>::n_dofs() const
{
  return dof_handler->n_dofs();
}

template <int dim>
void Heat<dim>::output_step(const double t,
                            const VEC &solution,
                            const VEC &solution_dot,
                            const unsigned int step_number,
                            const double /* h */ )
{
  computing_timer.enter_section ("Postprocessing");
  dirichlet_bcs.set_time(t);
  forcing_term.set_time(t);
  exact_solution.set_time(t);
  constraints.clear();
  DoFTools::make_hanging_node_constraints (*dof_handler,
                                           constraints);

  dirichlet_bcs.interpolate_boundary_values(*dof_handler, constraints);

  constraints.close ();
  VEC tmp(solution);
  constraints.distribute(tmp);
  distributed_solution = tmp;
  distributed_solution_dot = solution_dot;

  std::stringstream suffix;
  suffix << "." << current_cycle << "." << step_number;
  data_out.prepare_data_output( *dof_handler,
                                suffix.str());
  data_out.add_data_vector (distributed_solution, "u");
  std::vector<std::string> sol_dot_names =
    Utilities::split_string_list( "u");
  for (auto &name : sol_dot_names)
    {
      name += "_dot";
    }
  data_out.add_data_vector (distributed_solution_dot, print(sol_dot_names,","));

  data_out.write_data_and_clear("",*mapping);
  eh.error_from_exact(*mapping, *dof_handler, distributed_solution, exact_solution);

  computing_timer.exit_section ();
}

template <int dim>
bool Heat<dim>::solver_should_restart (const double t,
                                       const unsigned int step_number,
                                       const double h,
                                       VEC &solution,
                                       VEC &solution_dot)
{

  if (use_space_adaptivity)
    {
      int check = 0;
      double max_kelly=0;
      double mpi_max_kelly=0;

      computing_timer.enter_section ("   Compute error estimator");
      VEC tmp_c(solution);
      constraints.distribute(tmp_c);
      distributed_solution = tmp_c;

      Vector<float> estimated_error_per_cell (triangulation->n_active_cells());
      KellyErrorEstimator<dim>::estimate (*dof_handler,
                                          QGauss<dim-1>(fe->degree+1),
                                          typename FunctionMap<dim>::type(),
                                          distributed_solution,
                                          estimated_error_per_cell,
                                          ComponentMask(),
                                          0,
                                          0,
                                          triangulation->locally_owned_subdomain());


      max_kelly = estimated_error_per_cell.linfty_norm();
      max_kelly = Utilities::MPI::max(max_kelly, comm);

      if (max_kelly > kelly_threshold)

        {
          pcout << "  ################ restart ######### \n"
                << "max_kelly > threshold\n"
                << max_kelly  << " >  " << kelly_threshold
                << std::endl
                << "######################################\n";

          pgr.mark_cells(estimated_error_per_cell, *triangulation);

          parallel::distributed::SolutionTransfer<dim,VEC> sol_tr(*dof_handler);
          parallel::distributed::SolutionTransfer<dim,VEC> sol_dot_tr(*dof_handler);

          VEC sol (distributed_solution);
          VEC sol_dot (distributed_solution_dot);
          sol = solution;
          sol_dot = solution_dot;

          triangulation->prepare_coarsening_and_refinement();
          sol_tr.prepare_for_coarsening_and_refinement (sol);
          sol_dot_tr.prepare_for_coarsening_and_refinement(sol_dot);

          if (adaptive_refinement)
            triangulation->execute_coarsening_and_refinement ();
          else
            triangulation->refine_global (1);

          setup_dofs(false);

          VEC tmp (solution);
          VEC tmp_dot (solution_dot);

          sol_tr.interpolate (tmp);
          sol_dot_tr.interpolate (tmp_dot);

          solution = tmp;
          solution_dot = tmp_dot;
          constraints.distribute(solution);
          computing_timer.exit_section();
          MPI::COMM_WORLD.Barrier();
          return true;
        }
      else // if max_kelly > kelly_threshold
        {
          computing_timer.exit_section();
          return false;
        }

    }
  else // use space adaptivity

    return false;
}


template <int dim>
int Heat<dim>::setup_jacobian (const double t,
                               const VEC &src_yy,
                               const VEC &src_yp,
                               const VEC &,
                               const double alpha)
{
  computing_timer.enter_section ("   Setup Jacobian");
  assemble_jacobian_matrix(t, src_yy, src_yp, alpha);

//  TrilinosWrappers::PreconditionAMG::AdditionalData data;
//
//  preconditioner.initialize(jacobian_matrix, data);

  computing_timer.exit_section();

  return 0;
}

template <int dim>
int Heat<dim>::solve_jacobian_system (const double t,
                                      const VEC &y,
                                      const VEC &y_dot,
                                      const VEC &,
                                      const double alpha,
                                      const VEC &src,
                                      VEC &dst) const
{
  computing_timer.enter_section ("   Solve system");
  set_constrained_dofs_to_zero(dst);

  dst = Ainv * src;

  set_constrained_dofs_to_zero(dst);

  computing_timer.exit_section();
  return 0;
}

template <int dim>
VEC &Heat<dim>::differential_components() const
{
  static VEC diff_comps;
  IndexSet is = dof_handler->locally_owned_dofs();
  diff_comps.reinit(is, comm);
  diff_comps = 1;
  set_constrained_dofs_to_zero(diff_comps);
  return diff_comps;
}

template <int dim>
void Heat<dim>::set_constrained_dofs_to_zero(VEC &v) const
{
  for (unsigned int i=0; i<partitioning.n_elements(); ++i)
    {
      auto j = partitioning.nth_index_in_set(i);
      if (constraints.is_constrained(j))
        v[j] = 0;
    }
}

template <int dim>
void Heat<dim>::run ()
{

  for (current_cycle=0; current_cycle<n_cycles; ++current_cycle)
    {
      if (current_cycle == 0)
        {
          make_grid_fe();
          setup_dofs(true);
        }
      //  else
      //    refine_mesh();

      constraints.distribute(solution);

      dae.start_ode(solution, solution_dot, max_time_iterations);
      eh.error_from_exact(*mapping, *dof_handler, distributed_solution, exact_solution);
    }

  eh.output_table(pcout);

  // std::ofstream f("errors.txt");
  computing_timer.print_summary();
  timer_outfile.close();
  // f.close();
}

template class Heat<2>;

#endif
