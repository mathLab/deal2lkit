//-----------------------------------------------------------
//
//    Copyright (C) 2015 - 2016 by the deal2lkit authors
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


#include <deal2lkit/sundials_interface.h>
#include <deal2lkit/imex_stepper.h>
#ifdef D2K_WITH_SUNDIALS

#include <deal.II/base/utilities.h>
#include <deal.II/lac/block_vector.h>
#ifdef DEAL_II_WITH_TRILINOS
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#endif
#ifdef DEAL_II_WITH_PETSC
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#endif
#include <deal.II/base/utilities.h>

#include <iostream>
#include <iomanip>

#include <math.h>

#ifdef DEAL_II_WITH_MPI
#include <nvector/nvector_parallel.h>
#endif

using namespace dealii;


D2K_NAMESPACE_OPEN


template <typename VEC>
IMEXStepper<VEC>::IMEXStepper(SundialsInterface<VEC> &interface,
                              const double &step_size,
                              const double &initial_time,
                              const double &final_time) :
  ParameterAcceptor("IMEX Parameters"),
  interface(interface),
#ifdef DEAL_II_WITH_MPI
  kinsol("KINSOL for IMEX",interface.get_comm()),
#else
  kinsol("KINSOL for IMEX"),
#endif
  step_size(step_size),
  initial_time(initial_time),
  final_time(final_time),
  pcout(std::cout, Utilities::MPI::this_mpi_process(interface.get_comm())==0)
{
  abs_tol = 1e-6;
  rel_tol = 1e-8;
  output_period = 1;
  newton_alpha = 1.0;
  max_outer_non_linear_iterations = 5;
  max_inner_non_linear_iterations = 3;
  verbose = false;
  n_max_backtracking = 5;
  method = "fixed alpha";
}

template <typename VEC>
double IMEXStepper<VEC>::get_alpha() const
{
  double alpha;
  if (initial_time == final_time)
    alpha = 0.0;
  else
    alpha = 1./step_size;

  return alpha;

}

template <typename VEC>
void IMEXStepper<VEC>::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &step_size,
                "Step size", std::to_string(step_size), Patterns::Double());

  add_parameter(prm, &abs_tol,
                "Absolute error tolerance", std::to_string(abs_tol),
                Patterns::Double());

  add_parameter(prm, &rel_tol,
                "Relative error tolerance", std::to_string(rel_tol),
                Patterns::Double());

  add_parameter(prm, &initial_time,
                "Initial time", std::to_string(initial_time),
                Patterns::Double());

  add_parameter(prm, &final_time,
                "Final time", std::to_string(final_time),
                Patterns::Double());

  add_parameter(prm, &output_period,
                "Intervals between outputs", std::to_string(output_period),
                Patterns::Integer());

  add_parameter(prm, &max_outer_non_linear_iterations,
                "Maximum number of outer nonlinear iterations", std::to_string(max_outer_non_linear_iterations),
                Patterns::Integer(),
                "At each outer iteration the Jacobian is updated if it is set that the \n"
                "Jacobian is continuously updated and a cycle of inner iterations is \n"
                "perfomed.");

  add_parameter(prm, &max_inner_non_linear_iterations,
                "Maximum number of inner nonlinear iterations", std::to_string(max_inner_non_linear_iterations),
                Patterns::Integer(),
                "At each inner iteration the Jacobian is NOT updated.");

  add_parameter(prm, &newton_alpha,
                "Newton relaxation parameter", std::to_string(newton_alpha),
                Patterns::Double());

  add_parameter(prm, &update_jacobian_continuously,
                "Update continuously Jacobian", "true",
                Patterns::Bool());

  add_parameter(prm, &n_max_backtracking,
                "Number of elements in backtracking sequence", std::to_string(n_max_backtracking),
                Patterns::Integer(1),
                "In the line seach method with backtracking the following alphas are\n"
                "tested: 1, 1/2, 1/4,..., 2^-i. This parameter sets the maximum i.");
  add_parameter(prm, &method,
                "Method used", "fixed_alpha",
                Patterns::Selection("fixed_alpha|LS_backtracking"),
                "Fixed alpha means that the parsed alpha is used in each Newton iteration\n"
                "LS_backtracking is the line search with backtracking method.");

  add_parameter(prm, &verbose,
                "Print useful informations", "false",
                Patterns::Bool());

  add_parameter(prm, &use_kinsol,
                "Use the KINSOL solver", "true",
                Patterns::Bool());

}

template <typename VEC>
void IMEXStepper<VEC>::compute_y_dot(const VEC &y, const VEC &prev, const double alpha, VEC &y_dot)
{
  y_dot = y;
  y_dot -= prev;
  y_dot *= alpha;
}

template <typename VEC>
unsigned int IMEXStepper<VEC>::start_ode(VEC &solution, VEC &solution_dot)
{
  AssertDimension(solution.size(), interface.n_dofs());

  unsigned int step_number = 0;

  auto previous_solution = interface.create_new_vector();
  auto residual = interface.create_new_vector();
  auto rhs = interface.create_new_vector();

  double t = initial_time;
  double alpha;

  // check if it is a stationary problem
  if (initial_time == final_time)
    alpha = 0.0;
  else
    alpha = 1./step_size;

  compute_previous_solution(solution,solution_dot,alpha, *previous_solution);

  interface.output_step( 0, solution, solution_dot, 0, step_size);

  std::function<shared_ptr<VEC>()> my_new_vector = [this] ()
  {
    return this->interface.create_new_vector();
  };

  std::function<int(const VEC &, VEC &)> my_residual = [&] (const VEC &y, VEC &res)
  {
    compute_y_dot(y,*previous_solution,alpha,solution_dot);
    int ret = this->interface.residual(t,y,solution_dot,res);
    *residual = res;
    return ret;
  };

  std::function<int(const VEC &)> my_jac = [&] (const VEC &y)
  {
    compute_y_dot(y,*previous_solution,alpha,solution_dot);
    return this->interface.setup_jacobian(t,y,solution_dot,*residual,alpha);
  };

  std::function<int(const VEC &, VEC &)> my_solve = [&] (const VEC &res, VEC &dst)
  {
    *rhs = *residual;
    *rhs *= -1.0;
    return this->interface.solve_jacobian_system(t,solution,solution_dot,res,alpha,*rhs,dst);
  };

  std::function<int(const VEC &, VEC &)> my_jacobian_vmult = [&] (const VEC &v, VEC &dst )
  {
    return this->interface.jacobian_vmult(v,dst);
  };

  kinsol.create_new_vector = my_new_vector;
  kinsol.residual = my_residual;
  kinsol.setup_jacobian = my_jac;
  kinsol.solve_linear_system = my_solve;
  kinsol.jacobian_vmult = my_jacobian_vmult;

  // Initialization of the state of the boolean variable
  // responsible to keep track of the requirement that the
  // system's Jacobian be updated.
  bool update_Jacobian = true;

  // silence a warning when using kinsol
  (void) update_Jacobian;


  // store initial conditions
  interface.output_step(t, solution, solution_dot,  step_number, step_size);

  bool restart=false;

  auto L2 = interface.create_new_vector();
  if (use_kinsol)
    {
      // call kinsol initialization. this is mandatory if I am doing multiple cycle in pi-DoMUS
      kinsol.initialize_solver(solution);
      interface.get_lumped_mass_matrix(*L2);
      kinsol.set_scaling_vectors(*L2, *L2);
    }

  // The overall cycle over time begins here.
  while (t<=final_time+1e-15)
    {
      pcout << "Solving for t = " << t << std::endl;

      if (use_kinsol)
        {
          kinsol.solve(solution);
          compute_y_dot(solution,*previous_solution,alpha,solution_dot);
        }
      else
        do_newton(t,alpha,update_Jacobian,*previous_solution,solution,solution_dot);

      restart = interface.solver_should_restart(t,step_number,step_size,solution,solution_dot);

      while (restart)
        {

          previous_solution = interface.create_new_vector();
          residual = interface.create_new_vector();
          rhs = interface.create_new_vector();
          L2 = interface.create_new_vector();

          compute_previous_solution(solution,solution_dot,alpha,*previous_solution);

          if (use_kinsol)
            {


              kinsol.initialize_solver(solution);
              interface.get_lumped_mass_matrix(*L2);
              kinsol.set_scaling_vectors(*L2, *L2);
              kinsol.solve(solution);
              compute_y_dot(solution,*previous_solution,alpha,solution_dot);
            }
          else
            {
              do_newton(t,alpha,update_Jacobian,*previous_solution,solution,solution_dot);
            }

          restart = interface.solver_should_restart(t,step_number,step_size,solution,solution_dot);

        }


      step_number += 1;

      if ((step_number % output_period) == 0)
        interface.output_step(t, solution, solution_dot,  step_number, step_size);

      t += step_size;



      *previous_solution = solution;
      update_Jacobian = update_jacobian_continuously;

    } // End of the cycle over time.
  return 0;
}



template <typename VEC>
double IMEXStepper<VEC>::
line_search_with_backtracking(const VEC &update,
                              const VEC &previous_solution,
                              const double &alpha,
                              const double &t,
                              VEC &solution,
                              VEC &solution_dot,
                              VEC &residual)
{
  auto first_trial = interface.create_new_vector();
  auto first_residual = interface.create_new_vector();

  *first_trial = solution;
  double n_alpha = 1.0;

  first_trial->sadd(1.0, n_alpha, update);

  solution_dot = *first_trial;
  solution_dot -= previous_solution;
  solution_dot *= alpha;

  interface.residual(t, *first_trial, solution_dot, *first_residual);

  double first_res_norm = interface.vector_norm(*first_residual);

  auto second_trial = interface.create_new_vector();
  auto second_residual = interface.create_new_vector();

  for (unsigned int i=1; i<=n_max_backtracking; ++i)
    {
      *second_trial = solution;

      n_alpha = 1.0/(std::pow(2.0,i));

      second_trial->sadd(1.0, n_alpha, update);

      solution_dot = *second_trial;
      solution_dot -= previous_solution;
      solution_dot *= alpha;

      interface.residual(t, *second_trial, solution_dot, *second_residual);
      double second_res_norm = interface.vector_norm(*second_residual);
      if (first_res_norm < second_res_norm)
        {
          solution = *first_trial;
          residual = *first_residual;

          solution_dot = *first_trial;
          solution_dot -= previous_solution;
          solution_dot *= alpha;
          n_alpha = 1.0/(std::pow(2.0,i-1));
          return n_alpha;

        }
      else if (i<(n_max_backtracking))
        {
          *first_trial = *second_trial;
          *first_residual = *second_residual;
          first_res_norm = second_res_norm;
        }
      else
        {
          solution = *second_trial;
          residual = *second_residual;
          /* solution dot is ok */
          return n_alpha;
        }
    }
  return n_alpha;
}


template <typename VEC>
void IMEXStepper<VEC>::
do_newton (const double t,
           const double alpha,
           const bool update_Jacobian,
           const VEC &previous_solution,
           VEC &solution,
           VEC &solution_dot)
{
  auto solution_update = interface.create_new_vector();
  auto residual = interface.create_new_vector();
  auto rhs = interface.create_new_vector();



  // Initialization of two counters for the monitoring of
  // progress of the nonlinear solver.
  unsigned int inner_iter = 0;
  unsigned int outer_iter = 0;
  unsigned int nonlin_iter = 0;
  interface.residual(t, solution, solution_dot, *residual);
  double res_norm = 0.0;
  double solution_norm = 0.0;

  if (abs_tol>0.0||rel_tol>0.0)
    res_norm = interface.vector_norm(*residual);
  // if (rel_tol>0.0)
  //   solution_norm = interface.vector_norm(solution);

  // The nonlinear solver iteration cycle begins here.
  // using a do while approach, we ensure that the system
  // is solved at least once
  do
    {
      outer_iter += 1;
      if (update_Jacobian == true)
        {
          interface.setup_jacobian(t, solution, solution_dot,
                                   *residual, alpha);
        }

      inner_iter = 0;
      while (inner_iter < max_inner_non_linear_iterations &&
             res_norm > abs_tol &&
             res_norm > rel_tol*solution_norm)
        {


          inner_iter += 1;

          *rhs = *residual;
          *rhs *= -1.0;

          interface.solve_jacobian_system(t, solution, solution_dot,
                                          *residual, alpha,
                                          *rhs, *solution_update);


          if (method == "LS_backtracking")
            {
              newton_alpha = line_search_with_backtracking(*solution_update,
                                                           previous_solution,
                                                           alpha,
                                                           t,
                                                           solution,
                                                           solution_dot,
                                                           *residual);
            }
          else if (method == "fixed_alpha")
            {
              solution.sadd(1.0,
                            newton_alpha, *solution_update);

            }

          compute_y_dot(solution,previous_solution,alpha, solution_dot);

          res_norm = interface.vector_norm(*solution_update);

          if (rel_tol>0.0)
            {
              solution_norm = interface.vector_norm(solution);

              if (verbose)
                {
                  pcout << std::endl
                        << "   "
                        << " iteration "
                        << nonlin_iter + inner_iter
                        << ":\n"
                        << std::setw(19) << std::scientific << res_norm
                        << "   update norm\n"
                        << std::setw(19) << std::scientific << solution_norm
                        << "   solution norm\n"
                        << std::setw(19) << newton_alpha
                        << "   newton alpha\n\n"
                        << std::endl;
                }
            }
          else if (verbose)
            {
              pcout << std::endl
                    << "   "
                    << " iteration "
                    << nonlin_iter + inner_iter
                    << ":\n"
                    << std::setw(19) << std::scientific << res_norm
                    << "   update norm\n"
                    << std::setw(19) << newton_alpha
                    << "   newton alpha\n\n"
                    << std::endl;
            }

          interface.residual(t,solution,solution_dot,*residual);
        }

      nonlin_iter += inner_iter;

      if (std::fabs(res_norm) < abs_tol ||
          std::fabs(res_norm) < rel_tol*solution_norm)
        {
          pcout << std::endl
                << "   "
                << std::setw(19) << std::scientific << res_norm
                << " (converged in "
                << nonlin_iter
                << " iterations)\n\n"
                << std::endl;
          break; // Break of the while cycle
        }
      else if (outer_iter == max_outer_non_linear_iterations)
        {
          pcout << std::endl
                << "   "
                << std::setw(19) << std::scientific << res_norm
                << " (not converged in "
                << std::setw(3) << nonlin_iter
                << " iterations)\n\n"
                << std::endl;
          AssertThrow(false,
                      ExcMessage ("No convergence in nonlinear solver"));
        }

    } // The nonlinear solver iteration cycle ends here.
  while (outer_iter < max_outer_non_linear_iterations &&
         res_norm > abs_tol &&
         res_norm > rel_tol*solution_norm);
}

template <typename VEC>
void IMEXStepper<VEC>::
compute_previous_solution(const VEC &sol,
                          const VEC &sol_dot,
                          const double &alpha,
                          VEC &prev)
{
  if (alpha > 0.0)
    {
      prev = sol_dot;
      prev /= (-1.0*alpha);
      prev += sol;
    }
  else
    {
      prev = sol;
      (void)sol_dot;
      (void)alpha;
    }
}

D2K_NAMESPACE_CLOSE

template class deal2lkit::IMEXStepper<BlockVector<double> >;

#ifdef DEAL_II_WITH_MPI

#ifdef DEAL_II_WITH_TRILINOS
template class deal2lkit::IMEXStepper<TrilinosWrappers::MPI::Vector>;
template class deal2lkit::IMEXStepper<TrilinosWrappers::MPI::BlockVector>;
#endif

#ifdef DEAL_II_WITH_PETSC
template class deal2lkit::IMEXStepper<PETScWrappers::MPI::Vector>;
template class deal2lkit::IMEXStepper<PETScWrappers::MPI::BlockVector>;
#endif

#endif

#endif
