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
#include <deal.II/base/utilities.h>

#include <iostream>
#include <iomanip>

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
  step_size(step_size),
  initial_time(initial_time),
  final_time(final_time)
{
  abs_tol = 1e-6;
  rel_tol = 1e-8;
  output_period = 1;
  newton_alpha = 1.0;
  max_outer_non_linear_iterations = 5;
  max_inner_non_linear_iterations = 3;
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
}


template <typename VEC>
unsigned int IMEXStepper<VEC>::start_ode(VEC &solution)
{
  AssertDimension(solution.size(), interface.n_dofs());

  unsigned int step_number = 0;


  auto previous_solution = interface.create_new_vector();
  auto solution_dot = interface.create_new_vector();
  auto solution_update = interface.create_new_vector();
  auto residual = interface.create_new_vector();
  auto rhs = interface.create_new_vector();

  *previous_solution = solution;

  double t = initial_time;
  double alpha;

  // check if it is a stationary problem
  if (initial_time == final_time)
    alpha = 0.0;
  else
    alpha = 1./step_size;

  interface.output_step( 0, solution, *solution_dot, 0, step_size);

  // Initialization of the state of the boolean variable
  // responsible to keep track of the requirement that the
  // system's Jacobian be updated.
  bool update_Jacobian = true;

  // The overall cycle over time begins here.
  for (; t<=final_time; t+= step_size, ++step_number)
    {
      std::cout << "Time = " << t << std::endl;
      // Implicit Euler scheme.
      *solution_dot = solution;
      *solution_dot -= *previous_solution;
      *solution_dot *= alpha;

      // Initialization of two counters for the monitoring of
      // progress of the nonlinear solver.
      unsigned int inner_iter = 0;
      unsigned int outer_iter = 0;
      unsigned int nonlin_iter = 0;
      interface.residual(t, solution, *solution_dot, *residual);
      double res_norm = 0.0;
      double solution_norm = 0.0;

      if (abs_tol>0.0||rel_tol>0.0)
        res_norm = interface.vector_norm(*residual);
      if (rel_tol>0.0)
        solution_norm = interface.vector_norm(solution);

      // The nonlinear solver iteration cycle begins here.
      while (outer_iter < max_outer_non_linear_iterations &&
             res_norm > abs_tol &&
             res_norm > rel_tol*solution_norm)
        {
          outer_iter += 1;
          if (update_Jacobian == true)
            {
              interface.setup_jacobian(t, solution, *solution_dot,
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

              interface.solve_jacobian_system(t, solution, *solution_dot,
                                              *residual, alpha,
                                              *rhs, *solution_update);
              solution.sadd(1.0,
                            newton_alpha, *solution_update);

              // Implicit Euler scheme.
              *solution_dot = solution;
              *solution_dot -= *previous_solution;
              *solution_dot *= alpha;

              interface.residual(t, solution, *solution_dot, *residual);

              if (abs_tol>0.0||rel_tol>0.0)
                res_norm = interface.vector_norm(*solution_update);
              if (rel_tol>0.0)
                solution_norm = interface.vector_norm(solution);

            }

          nonlin_iter += inner_iter;

          if (std::fabs(res_norm) < abs_tol)
            {
              std::printf("   %-16.3e (converged in %d iterations)\n\n", res_norm, nonlin_iter);
              break; // Break of the while cycle ... after this a time advancement happens.
            }
          else if (outer_iter == max_outer_non_linear_iterations)
            {
              std::printf("   %-16.3e (not converged in %d iterations)\n\n", res_norm, nonlin_iter);
              AssertThrow(false,
                          ExcMessage ("No convergence in nonlinear solver"));
            }

        } // The nonlinear solver iteration cycle ends here.

      *previous_solution = solution;

      if ((step_number % output_period) == 0)
        interface.output_step(t, solution, *solution_dot,  step_number, step_size);

      update_Jacobian = update_jacobian_continuously;

    } // End of the cycle over time.
  return 0;
}



D2K_NAMESPACE_CLOSE

template class deal2lkit::IMEXStepper<BlockVector<double> >;

#ifdef DEAL_II_WITH_MPI

#ifdef DEAL_II_WITH_TRILINOS
template class deal2lkit::IMEXStepper<TrilinosWrappers::MPI::Vector>;
template class deal2lkit::IMEXStepper<TrilinosWrappers::MPI::BlockVector>;
#endif

#endif

#endif
