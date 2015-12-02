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
}

template <typename VEC>
void IMEXStepper<VEC>::declare_parameters(ParameterHandler &prm)
{
add_parameter(prm, &step_size,
"Initial step size", "1e-4", Patterns::Double());

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

  add_parameter(prm, &max_non_linear_iterations,
                "Maximum number of nonlinear iterations", std::to_string(max_non_linear_iterations),
                Patterns::Integer());

add_parameter(prm, &newton_alpha,
              "Newton relaxation parameter", std::to_string(newton_alpha),
              Patterns::Double());
}


template <typename VEC>
unsigned int IMEXStepper<VEC>::start_ode(VEC &solution)
{
  AssertDimension(solution.size(), interface.n_dofs());

  unsigned int step_number = 0;

  int status;

    // The solution is stored in
    // solution. Here we take only a
    // view of it.

  auto previous_solution = interface.create_new_vector();
  auto solution_dot = interface.create_new_vector();
  auto solution_update = interface.create_new_vector();
  auto residual = interface.create_new_vector();

  *previous_solution = solution;

  double t = initial_time;
  double alpha = 1./step_size;

  interface.output_step( 0, solution, solution_dot, 0, step_size);

    // Initialization of the state of the boolean variable
    // responsible to keep track of the requirement that the
    // system's Jacobian be updated.
  bool update_Jacobian = true;

    // The overall cycle over time begins here.
  for(; t<=final_time; t+= step_size, ++step_number)
  {
      // Initialization of two counters for the monitoring of
      // progress of the nonlinear solver.
    unsigned int nonlin_iter = 0;
    unsigned int outer_nonlin_iter = 0;

      // The nonlinear solver iteration cycle begins here.
    while (true)
	 {
       // Implicit Euler scheme.
     solution_dot->sadd(0,
                        alpha, *solution,
                        -alpha, *previous_solution);

     xi = current_xi;

     if(update_Jacobian == true)
     {
       setup_jacobian(t, *solution, *solution_dot,
                      *residual, alpha);
     }

     interface.residual(t, *solution, *solution_dot, *residual);

     *residual *= -1.0;

     solve_jacobian_system(t, *solution, *solution_dot,
                           *residual, alpha,
                           *residual, *solution_update);

     update_Jacobian = update_jacobian_continuously;
     }
	    else
        assemble_residual (current_res, JF, xit, xi,
                           par.dt, t);

     const double res_norm = current_res.l2_norm();

     if (std::fabs(res_norm) < 1e-10)
     {
       std::printf("   %-16.3e (converged in %d iterations)\n\n", res_norm, nonlin_iter);
       break; // Break of the while cycle ... after this a time advancement happens.
     }
     else
     {
       current_res *= -1;
       static Vector<double> tmp(current_res.size());
       tmp = current_res;
       JF_inv.solve(tmp);
       newton_update = tmp;
       constraints_f.distribute(newton_update.block(0));
       current_xi.add(1., newton_update);
       double avg_pressure = (current_xi.block(0)*pressure_average);
       current_xi.block(0).add(-avg_pressure, unit_pressure);
       if(res_norm > 1e-3) update_Jacobian = true;
     }

     ++nonlin_iter;
     if(nonlin_iter == 10)
     {
       update_Jacobian = true;
       nonlin_iter = 0;
       outer_nonlin_iter++;
       std::printf("   %-16.3e (not converged in 10 iterations. Step %d)\n\n", res_norm, outer_nonlin_iter);
     }

     AssertThrow (outer_nonlin_iter <= 3,
                  ExcMessage ("No convergence in nonlinear solver"));

   } // The nonlinear solver iteration cycle ends here.
  previous_xi = current_xi;

  output_results(current_xi, time_step);
  update_Jacobian = true;
  
} // End of the cycle over time.
}
}
status = IDASolve(ida_mem, next_time, &t, yy, yp, IDA_NORMAL);

status = IDAGetLastStep(ida_mem, &h);
AssertThrow(status == 0, ExcMessage("Error in IDA Solver"));

copy(solution, yy);
copy(solution_dot, yp);

// Check the solution
bool reset = solver.solver_should_restart(t, step_number, h, solution, solution_dot);


while ( reset == true )
{
// double frac = 0;
int k = 0;
IDAGetLastOrder(ida_mem, &k);
// frac = std::pow((double)k,2.);
reset_ode(t, solution, solution_dot,
h/2.0, max_steps, false);
reset = solver.solver_should_restart(t, step_number, h, solution, solution_dot);
}

step_number++;

solver.output_step(t, solution, solution_dot,  step_number, h);



}

std::cout << std::endl;
// Free the vectors which are no longer used.
#ifdef DEAL_II_WITH_MPI
N_VDestroy_Parallel(yy);
N_VDestroy_Parallel(yp);
N_VDestroy_Parallel(abs_tolls);
N_VDestroy_Parallel(diff_id);
#else
N_VDestroy_Serial(yy);
N_VDestroy_Serial(yp);
N_VDestroy_Serial(abs_tolls);
N_VDestroy_Serial(diff_id);
#endif

return step_number;
}

template <typename VEC>
void IMEXStepper<VEC>::reset_ode(double current_time,
VEC &solution,
VEC &solution_dot,
double current_time_step,
unsigned int max_steps,
bool first_step)
{
if (ida_mem)
IDAFree(&ida_mem);

ida_mem = IDACreate();


// Free the vectors which are no longer used.
if (yy)
{
#ifdef DEAL_II_WITH_MPI
N_VDestroy_Parallel(yy);
N_VDestroy_Parallel(yp);
N_VDestroy_Parallel(abs_tolls);
N_VDestroy_Parallel(diff_id);
#else
N_VDestroy_Serial(yy);
N_VDestroy_Serial(yp);
N_VDestroy_Serial(abs_tolls);
N_VDestroy_Serial(diff_id);
#endif
}

int status;
Assert(solution.size() == solver.n_dofs(),
ExcDimensionMismatch(solution.size(), solver.n_dofs()));

Assert(solution_dot.size() == solver.n_dofs(),
ExcDimensionMismatch(solution_dot.size(), solver.n_dofs()));


#ifdef DEAL_II_WITH_MPI
IndexSet is = solution.locally_owned_elements();

yy        = N_VNew_Parallel(solver.get_comm(), is.n_elements(), solver.n_dofs());
yp        = N_VNew_Parallel(solver.get_comm(), is.n_elements(), solver.n_dofs());
diff_id   = N_VNew_Parallel(solver.get_comm(), is.n_elements(), solver.n_dofs());
abs_tolls = N_VNew_Parallel(solver.get_comm(), is.n_elements(), solver.n_dofs());
#else
yy        = N_VNew_Serial(solver.n_dofs());
yp        = N_VNew_Serial(solver.n_dofs());
diff_id   = N_VNew_Serial(solver.n_dofs());
abs_tolls = N_VNew_Serial(solver.n_dofs());
#endif

copy(yy, solution);
copy(yp, solution_dot);
copy(diff_id, solver.differential_components());

status = IDAInit(ida_mem, t_dae_residual<VEC>, current_time, yy, yp);

if (use_local_tolerances)
{
VEC &tolerances = solver.get_local_tolerances();
VEC abs_tolerances(tolerances);
abs_tolerances /= tolerances.linfty_norm();
abs_tolerances *= abs_tol;
copy(abs_tolls, abs_tolerances);
status += IDASVtolerances(ida_mem, rel_tol, abs_tolls);
}
else
{
status += IDASStolerances(ida_mem, rel_tol, abs_tol);
}

status += IDASetInitStep(ida_mem, current_time_step);
status += IDASetUserData(ida_mem, (void *) &solver);

status += IDASetId(ida_mem, diff_id);
status += IDASetSuppressAlg(ida_mem, ignore_algebraic_terms_for_errors);

status += IDASetMaxNumSteps(ida_mem, max_steps);
status += IDASetStopTime(ida_mem, final_time);

status += IDASetMaxNonlinIters(ida_mem, max_non_linear_iterations);

// Initialize solver
IDAMem IDA_mem;
IDA_mem = (IDAMem) ida_mem;

IDA_mem->ida_lsetup = t_dae_lsetup<VEC>;
IDA_mem->ida_lsolve = t_dae_solve<VEC>;
IDA_mem->ida_setupNonNull = true;


status += IDASetMaxOrd(ida_mem, max_order);

AssertThrow(status == 0, ExcMessage("Error initializing IDA."));

std::string type;
if (first_step)
type = ic_type;
else
type = reset_type;

if (type == "use_y_dot")
{
// (re)initialization of the vectors
IDACalcIC(ida_mem, IDA_Y_INIT, current_time+current_time_step);
IDAGetConsistentIC(ida_mem, yy, yp);

copy(solution, yy);
copy(solution_dot, yp);
}
else if (type == "use_y_diff")
{
IDACalcIC(ida_mem, IDA_YA_YDP_INIT, current_time+current_time_step);
IDAGetConsistentIC(ida_mem, yy, yp);

copy(solution, yy);
copy(solution_dot, yp);
}
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
