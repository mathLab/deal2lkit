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

#include "stokes_ida.h"

#ifdef D2K_WITH_SUNDIALS

#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_solver.h>

#include <nvector/nvector_parallel.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/numerics/solution_transfer.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_math.h>

#include <deal.II/lac/solver_gmres.h>
template <int dim>
Stokes<dim>::Stokes (const MPI_Comm &communicator)
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

  pgg("Domain"),
  pgr("Refinement"),
  fe_builder("Finite Element",
             "FESystem[FE_Q(2)^dim-FE_Q(1)]",
             "u,u,p"),

  exact_solution("Exact solution", 3),
  forcing_term("Forcing term", 3),
  initial_solution("Initial solution", 3),
  initial_solution_dot("Initial solution_dot", 3),
  dirichlet_bcs("Dirichlet BCs", 3, "u,u,p", "0=u"),
  dirichlet_dot("Dirichlet dot", 3, "u,u,p", "0=u"),

  data_out("Output Parameters", "vtu"),
  dae(*this)
{}

template <int dim>
void Stokes<dim>::declare_parameters (ParameterHandler &prm)
{
  add_parameter(  prm,
                  &initial_global_refinement,
                  "Initial global refinement",
                  "1",
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
                  &mu,
                  "mu",
                  "1.",
                  Patterns::Double(0.0));
}

template <int dim>
void Stokes<dim>::make_grid_fe()
{
  fe=SP(fe_builder());
  triangulation = SP(pgg.distributed(comm));

  static HyperShellBoundary<dim> boundary;
  triangulation->set_all_manifold_ids_on_boundary(5,99);
  triangulation->set_manifold (99, boundary);
  triangulation->refine_global (initial_global_refinement);
  dof_handler = SP(new DoFHandler<dim>(*triangulation));
}


template <int dim>
void Stokes<dim>::setup_dofs (const bool &first_run)
{
  computing_timer.enter_section("Setup dof systems");

  std::vector<unsigned int> sub_blocks = fe_builder.get_component_blocks();

  dof_handler->distribute_dofs (*fe);
  DoFRenumbering::component_wise (*dof_handler, sub_blocks);

  mapping = SP(new MappingQ<dim>(2));

  dofs_per_block.clear();
  dofs_per_block.resize(fe_builder.n_blocks());

  DoFTools::count_dofs_per_block (*dof_handler, dofs_per_block,
                                  sub_blocks);

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
        << " (" << print(dofs_per_block, "+") << ")"
        << std::endl;
  pcout.get_stream().imbue(s);


  partitioning.resize(0);
  relevant_partitioning.resize(0);

  IndexSet relevant_set;
  {
    global_partitioning = dof_handler->locally_owned_dofs();
    for (unsigned int i = 0; i < fe_builder.n_blocks(); ++i)
      partitioning.push_back(global_partitioning.get_view( std::accumulate(dofs_per_block.begin(), dofs_per_block.begin() + i, 0),
                                                           std::accumulate(dofs_per_block.begin(), dofs_per_block.begin() + i + 1, 0)));

    DoFTools::extract_locally_relevant_dofs (*dof_handler,
                                             relevant_set);

    for (unsigned int i = 0; i < fe_builder.n_blocks(); ++i)
      relevant_partitioning.push_back(relevant_set.get_view(std::accumulate(dofs_per_block.begin(), dofs_per_block.begin() + i, 0),
                                                            std::accumulate(dofs_per_block.begin(), dofs_per_block.begin() + i + 1, 0)));
  }
  constraints.clear ();
  constraints.reinit (relevant_set);

  DoFTools::make_hanging_node_constraints (*dof_handler,
                                           constraints);

  dirichlet_bcs.interpolate_boundary_values(*dof_handler, constraints);
  constraints.close ();

  constraints_dot.clear ();
  constraints_dot.reinit (relevant_set);

  DoFTools::make_hanging_node_constraints (*dof_handler,
                                           constraints_dot);

  dirichlet_dot.interpolate_boundary_values(*dof_handler, constraints_dot);
  constraints_dot.close ();

  jacobian_matrix.clear();
  jacobian_matrix_sp.reinit(partitioning,partitioning,relevant_partitioning,comm);


  DoFTools::make_sparsity_pattern (*dof_handler,
                                   fe_builder.get_coupling(),
                                   jacobian_matrix_sp,
                                   constraints,
                                   false,
                                   Utilities::MPI::this_mpi_process(comm));

  jacobian_matrix_sp.compress();

  jacobian_matrix.reinit(jacobian_matrix_sp);

  jacobian_preconditioner_matrix.clear();
  jacobian_preconditioner_matrix_sp.reinit(partitioning,partitioning,relevant_partitioning,comm);


  DoFTools::make_sparsity_pattern (*dof_handler,
                                   fe_builder.get_preconditioner_coupling(),
                                   jacobian_preconditioner_matrix_sp,
                                   constraints,
                                   false,
                                   Utilities::MPI::this_mpi_process(comm));

  jacobian_preconditioner_matrix_sp.compress();

  jacobian_preconditioner_matrix.reinit(jacobian_preconditioner_matrix_sp);

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
void Stokes<dim>::update_constraints (const double &t)
{
  dirichlet_bcs.set_time(t);
  dirichlet_dot.set_time(t);
  exact_solution.set_time(t);
  constraints.clear();
  DoFTools::make_hanging_node_constraints (*dof_handler,
                                           constraints);

  dirichlet_bcs.interpolate_boundary_values(*dof_handler, constraints);

  constraints.close ();
  constraints_dot.clear();

  DoFTools::make_hanging_node_constraints (*dof_handler,
                                           constraints_dot);

  dirichlet_dot.interpolate_boundary_values(*dof_handler, constraints_dot);
  constraints_dot.close ();
}

template <int dim>
void Stokes<dim>::assemble_jacobian_matrix(const double t,
                                           const VEC &solution,
                                           const VEC &solution_dot,
                                           const double alpha)
{
  computing_timer.enter_section ("   Assemble jacobian matrix");
  jacobian_matrix = 0;
  jacobian_preconditioner_matrix = 0;

  update_constraints(t);

  VEC tmp(solution);
  VEC tmp_dot(solution_dot);
  constraints.distribute(tmp);
  constraints_dot.distribute(tmp_dot);
  distributed_solution = tmp;
  distributed_solution_dot = tmp_dot;

  const QGauss<dim>  quadrature_formula(fe->degree+1);

  FEValues<dim> fe_values (*mapping,*fe, quadrature_formula,
                           update_values    |  update_gradients |
                           update_quadrature_points |
                           update_JxW_values);

  const unsigned int   dofs_per_cell = fe->dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double>   cell_prec (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  std::vector<Tensor<1,dim> >          phi_u(dofs_per_cell);
  std::vector<SymmetricTensor<2,dim> > sym_grads_phi_u(dofs_per_cell);
  std::vector<Tensor<2,dim> >          grads_phi_u(dofs_per_cell);
  std::vector<double>                  div_phi_u(dofs_per_cell);
  std::vector<double>                  phi_p(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler->begin_active(),
  endc = dof_handler->end();
  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned())
      {
        cell_matrix = 0;
        cell_prec = 0;

        fe_values.reinit (cell);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          {
            for (unsigned int k=0; k<dofs_per_cell; ++k)
              {
                phi_u[k]       = fe_values[velocities].value (k,q_point);
                sym_grads_phi_u[k] = fe_values[velocities].symmetric_gradient(k,q_point);
                grads_phi_u[k] = fe_values[velocities].gradient(k,q_point);
                div_phi_u[k]   = fe_values[velocities].divergence (k, q_point);
                phi_p[k]       = fe_values[pressure].value (k, q_point);
              }

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  {
                    cell_matrix(i,j) += (
                                          alpha*phi_u[i] *
                                          phi_u[j]
                                          +
                                          mu*sym_grads_phi_u[i] * sym_grads_phi_u[j]
                                          -
                                          (div_phi_u[i] * phi_p[j])
                                          -
                                          (phi_p[i] * div_phi_u[j])

                                        )*fe_values.JxW(q_point);

                    cell_prec(i,j) += (
                                        alpha*phi_u[i]*phi_u[j]
                                        +
                                        mu*sym_grads_phi_u[i]*sym_grads_phi_u[j]
                                        +
                                        (1./mu)*phi_p[i]*phi_p[j]

                                      )*fe_values.JxW(q_point);
                  }


              }
          }

        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (cell_matrix,
                                                local_dof_indices,
                                                jacobian_matrix);

        constraints.distribute_local_to_global (cell_prec,
                                                local_dof_indices,
                                                jacobian_preconditioner_matrix);
      }

  jacobian_matrix.compress (VectorOperation::add);
  jacobian_preconditioner_matrix.compress (VectorOperation::add);

  // ########## Operators setup
  std::vector<std::vector<bool> > constant_modes;
  DoFTools::extract_constant_modes (*dof_handler, dof_handler->get_fe().component_mask(velocities),
                                    constant_modes);

  Mp_preconditioner.reset  (new TrilinosWrappers::PreconditionJacobi());
  Amg_preconditioner.reset (new TrilinosWrappers::PreconditionAMG());

  TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data;
  Amg_data.constant_modes = constant_modes;
  Amg_data.elliptic = true;
  Amg_data.higher_order_elements = true;
  Amg_data.smoother_sweeps = 2;
  Amg_data.aggregation_threshold = 0.02;

  Mp_preconditioner->initialize  (jacobian_preconditioner_matrix.block(1,1));
  Amg_preconditioner->initialize (jacobian_matrix.block(0,0),
                                  Amg_data);


  // SYSTEM MATRIX:
  std::array<std::array<LinearOperator<TrilinosWrappers::MPI::Vector>, 2>, 2> S;
  for (unsigned int i = 0; i<2; ++i)
    for (unsigned int j = 0; j<2; ++j)
      S[i][j] = linear_operator< TrilinosWrappers::MPI::Vector >( jacobian_matrix.block(i,j) );

  S[1][1] = null_operator(S[1][1]); // p p

  jacobian_op = BlockLinearOperator<TrilinosWrappers::MPI::BlockVector>(S);

  // PRECONDITIONER

  std::array<std::array<LinearOperator<TrilinosWrappers::MPI::Vector>, 2>, 2> P;
  for (unsigned int i = 0; i<2; ++i)
    for (unsigned int j = 0; j<2; ++j)
      P[i][j] = linear_operator< TrilinosWrappers::MPI::Vector >( jacobian_preconditioner_matrix.block(i,j) );

  static ReductionControl solver_control_pre(5000, 1e-8);
  static SolverCG<TrilinosWrappers::MPI::Vector> solver_CG(solver_control_pre);

  for (unsigned int i = 0; i<2; ++i)
    for (unsigned int j = 0; j<2; ++j)
      if (i!=j)
        P[i][j] = null_operator< TrilinosWrappers::MPI::Vector >(P[i][j]);

  P[0][0] = inverse_operator< >(  S[0][0],
                                  solver_CG,
                                  *Amg_preconditioner);
  P[1][1] = inverse_operator< >(  jacobian_preconditioner_matrix.block(1,1),
                                  solver_CG,
                                  *Mp_preconditioner);

  jacobian_preconditioner_op = block_forward_substitution< >(
                                 BlockLinearOperator<TrilinosWrappers::MPI::BlockVector>(S),
                                 BlockLinearOperator<TrilinosWrappers::MPI::BlockVector>(P));

  computing_timer.exit_section();
}

template <int dim>
int Stokes<dim>::residual (const double t,
                           const VEC &solution,
                           const VEC &solution_dot,
                           VEC &dst)
{
  computing_timer.enter_section ("Residual");

  update_constraints(t);


  VEC tmp(solution);
  VEC tmp_dot(solution_dot);
  constraints.distribute(tmp);
  constraints_dot.distribute(tmp_dot);

  distributed_solution = tmp;
  distributed_solution_dot = tmp_dot;

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

  const FEValuesExtractors::Vector u (0);
  const FEValuesExtractors::Scalar p (dim);

  std::vector<Tensor<1,dim> >          phi_u(dofs_per_cell);
  std::vector<SymmetricTensor<2,dim> > grads_phi_u(dofs_per_cell);
  std::vector<double>                  div_phi_u(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler->begin_active(),
  endc = dof_handler->end();
  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned())
      {
        cell_rhs = 0;

        fe_values.reinit (cell);
        cell->get_dof_indices (local_dof_indices);

        std::vector<Vector<double> > rhs_values (n_q_points, Vector<double>(dim+1));
        forcing_term.vector_value_list (fe_values.get_quadrature_points(),
                                        rhs_values);


        std::vector<SymmetricTensor<2,dim> > grad_sols(n_q_points);
        fe_values[u].get_function_symmetric_gradients(distributed_solution,grad_sols);

        std::vector<Tensor<1,dim> > sols_dot (n_q_points);
        fe_values[u].get_function_values(distributed_solution_dot,sols_dot);

        std::vector<double> div_us (n_q_points);
        fe_values[u].get_function_divergences(distributed_solution,div_us);

        std::vector<double> ps (n_q_points);
        fe_values[p].get_function_values(distributed_solution,ps);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          {

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {

                cell_rhs(i) += (
                                 sols_dot[q_point]*
                                 fe_values[u].value(i,q_point)

                                 +
                                 mu*scalar_product(grad_sols[q_point],
                                                   fe_values[u].symmetric_gradient(i,q_point))

                                 -
                                 ps[q_point] *
                                 fe_values[u].divergence(i,q_point)

                                 -
                                 div_us[q_point]*
                                 fe_values[p].value(i,q_point)

                               )*fe_values.JxW(q_point);


                unsigned int comp_i = fe->system_to_component_index(i).first;
                if (comp_i<dim)
                  cell_rhs(i) -= (rhs_values[q_point](comp_i) *
                                  fe_values[u].value(i,q_point)[comp_i] *
                                  fe_values.JxW(q_point));
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
shared_ptr<VEC> Stokes<dim>::create_new_vector () const
{
  shared_ptr<VEC> ret = SP(new VEC(solution));
  return ret;
}

template <int dim>
unsigned int Stokes<dim>::n_dofs() const
{
  return dof_handler->n_dofs();
}

template <int dim>
void Stokes<dim>::output_step(const double t,
                              const VEC &solution,
                              const VEC &solution_dot,
                              const unsigned int step_number,
                              const double /* h */ )
{
  computing_timer.enter_section ("Postprocessing");

  update_constraints(t);

  VEC tmp(solution);
  VEC tmp_dot(solution_dot);

  constraints.distribute(tmp);
  constraints_dot.distribute(tmp_dot);
  distributed_solution = tmp;
  distributed_solution_dot = tmp_dot;

  std::stringstream suffix;
  suffix << "." << step_number;
  data_out.prepare_data_output( *dof_handler,
                                suffix.str());
  data_out.add_data_vector (distributed_solution, fe_builder.get_component_names());
  std::vector<std::string> sol_dot_names =
    Utilities::split_string_list( fe_builder.get_component_names());
  for (auto &name : sol_dot_names)
    {
      name += "_dot";
    }
  data_out.add_data_vector (distributed_solution_dot, print(sol_dot_names,","));

  data_out.write_data_and_clear(*mapping);

  computing_timer.exit_section ();
}

template <int dim>
bool Stokes<dim>::solver_should_restart (const double t,
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


      update_constraints(t);

      VEC tmp_c(solution);
      constraints.distribute(tmp_c);
      distributed_solution = tmp_c;
      std::vector<bool> mask(dim+1,true);
      mask[dim] = false;

      Vector<float> estimated_error_per_cell (triangulation->n_active_cells());
      KellyErrorEstimator<dim>::estimate (*dof_handler,
                                          QGauss<dim-1>(4),
                                          typename FunctionMap<dim>::type(),
                                          distributed_solution,
                                          estimated_error_per_cell,
                                          ComponentMask(mask),
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

          update_constraints(t);

          constraints.distribute(solution);
          constraints_dot.distribute(solution_dot);
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
int Stokes<dim>::setup_jacobian (const double t,
                                 const VEC &src_yy,
                                 const VEC &src_yp,
                                 const VEC &,
                                 const double alpha)
{
  computing_timer.enter_section ("   Setup Jacobian");

  assemble_jacobian_matrix(t, src_yy, src_yp, alpha);

  computing_timer.exit_section();

  return 0;
}

template <int dim>
int Stokes<dim>::solve_jacobian_system (const double t,
                                        const VEC &y,
                                        const VEC &y_dot,
                                        const VEC &,
                                        const double alpha,
                                        const VEC &src,
                                        VEC &dst) const
{
  computing_timer.enter_section ("   Solve system");

  dst = solution;
  set_constrained_dofs_to_zero(dst);


//  const double solver_tolerance = 1e-8*src.l2_norm();
//
//  PrimitiveVectorMemory<VEC> mem;
//  SolverControl solver_control (30, solver_tolerance);
//  SolverControl solver_control_refined (jacobian_matrix.m(), solver_tolerance);
//
//  SolverFGMRES<VEC>
//  solver(solver_control, mem,
//         typename SolverFGMRES<VEC>::AdditionalData(30, true));
//
//  SolverFGMRES<VEC>
//  solver_refined(solver_control_refined, mem,
//                 typename SolverFGMRES<VEC>::AdditionalData(50, true));
//
//  auto S_inv         = inverse_operator(jacobian_op, solver, jacobian_preconditioner_op);
//  auto S_inv_refined = inverse_operator(jacobian_op, solver_refined, jacobian_preconditioner_op);
//  unsigned int n_iterations = 0;
//  try
//    {
//      S_inv.vmult(dst, src);
//      n_iterations = solver_control.last_step();
//    }
//  catch ( SolverControl::NoConvergence )
//    {
//      try
//        {
//          S_inv_refined.vmult(dst, src);
//          n_iterations = (solver_control.last_step() +
//                          solver_control_refined.last_step());
//        }
//      catch ( SolverControl::NoConvergence )
//        {
//          computing_timer.exit_section();
//          return 1;
//        }
//
//    }
//  pcout << std::endl
//        << " iterations:                           " <<  n_iterations
//        << std::endl;
//
//  set_constrained_dofs_to_zero(dst);

  PrimitiveVectorMemory<TrilinosWrappers::MPI::BlockVector> mem;

  unsigned int n_iterations = 0;
  const double solver_tolerance = 1e-8 * src.l2_norm();
  SolverControl solver_control (30, solver_tolerance);

  try
    {
      const LinearSolvers::BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG,
            TrilinosWrappers::PreconditionJacobi>
            preconditioner (jacobian_matrix, jacobian_preconditioner_matrix,
                            *Mp_preconditioner, *Amg_preconditioner,
                            false);

      SolverFGMRES<TrilinosWrappers::MPI::BlockVector>
      solver(solver_control, mem,
             SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::
             AdditionalData(30, true));
      solver.solve(jacobian_matrix, dst, src,
                   preconditioner);

      n_iterations = solver_control.last_step();
    }

  catch (SolverControl::NoConvergence)
    {
      const LinearSolvers::BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG,
            TrilinosWrappers::PreconditionJacobi>
            preconditioner (jacobian_matrix, jacobian_preconditioner_matrix,
                            *Mp_preconditioner, *Amg_preconditioner,
                            true);

      SolverControl solver_control_refined (jacobian_matrix.m(), solver_tolerance);
      SolverFGMRES<TrilinosWrappers::MPI::BlockVector>
      solver(solver_control_refined, mem,
             SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::
             AdditionalData(50, true));
      solver.solve(jacobian_matrix, dst, src,
                   preconditioner);

      n_iterations = (solver_control.last_step() +
                      solver_control_refined.last_step());
    }


  set_constrained_dofs_to_zero(dst);
  pcout << n_iterations  << " iterations."
        << std::endl;
  computing_timer.exit_section();
  return 0;
}

template <int dim>
VEC &Stokes<dim>::differential_components() const
{
  static VEC diff_comps;
  diff_comps.reinit(solution);
  diff_comps.block(0) = 1; // velocity is differential
  diff_comps.block(1) = 0; // pressure is only algebraic

  set_constrained_dofs_to_zero(diff_comps);
  return diff_comps;
}

template <int dim>
void Stokes<dim>::set_constrained_dofs_to_zero(VEC &v) const
{
  for (unsigned int i=0; i<global_partitioning.n_elements(); ++i)
    {
      auto j = global_partitioning.nth_index_in_set(i);
      if (constraints.is_constrained(j))
        v[j] = 0;
    }
}

template <int dim>
void Stokes<dim>::run ()
{

  make_grid_fe();
  setup_dofs(true);

  constraints.distribute(solution);
  constraints_dot.distribute(solution_dot);

  dae.start_ode(solution, solution_dot, max_time_iterations);

  computing_timer.print_summary();
  timer_outfile.close();
}

template class Stokes<2>;

#endif
