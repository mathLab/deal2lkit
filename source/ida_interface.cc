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


#include <deal2lkit/ida_interface.h>
#include <deal2lkit/sundials_interface.h>
#include <deal2lkit/utilities.h>

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
#include <ida/ida_impl.h>


using namespace dealii;


D2K_NAMESPACE_OPEN


template<typename VEC>
int t_dae_residual(realtype tt, N_Vector yy, N_Vector yp,
                   N_Vector rr, void *user_data)
{
  SundialsInterface<VEC> &solver = *static_cast<SundialsInterface<VEC> *>(user_data);

  shared_ptr<VEC> src_yy = solver.create_new_vector();
  shared_ptr<VEC> src_yp = solver.create_new_vector();
  shared_ptr<VEC> residual = solver.create_new_vector();

  copy(*src_yy, yy);
  copy(*src_yp, yp);

  int err = solver.residual(tt, *src_yy, *src_yp, *residual);

  copy(rr, *residual);

  return err;
}



template<typename VEC>
int t_dae_lsetup(IDAMem IDA_mem,
                 N_Vector yy,
                 N_Vector yp,
                 N_Vector resp,
                 N_Vector tmp1,
                 N_Vector tmp2,
                 N_Vector tmp3)
{
  (void) tmp1;
  (void) tmp2;
  (void) tmp3;
  SundialsInterface<VEC> &solver = *static_cast<SundialsInterface<VEC> *>(IDA_mem->ida_user_data);

  shared_ptr<VEC> src_yy = solver.create_new_vector();
  shared_ptr<VEC> src_yp = solver.create_new_vector();
  shared_ptr<VEC> residual = solver.create_new_vector();

  copy(*src_yy, yy);
  copy(*src_yp, yp);
  copy(*residual, resp);

  int err = solver.setup_jacobian(IDA_mem->ida_tn,
                                  *src_yy,
                                  *src_yp,
                                  *residual,
                                  IDA_mem->ida_cj);
  return err;
}


template<typename VEC>
int t_dae_solve(IDAMem IDA_mem,
                N_Vector b,
                N_Vector weight,
                N_Vector yy,
                N_Vector yp,
                N_Vector resp)
{
  (void) weight;
  SundialsInterface<VEC> &solver = *static_cast<SundialsInterface<VEC> *>(IDA_mem->ida_user_data);

  shared_ptr<VEC> src_yy = solver.create_new_vector();
  shared_ptr<VEC> src_yp = solver.create_new_vector();
  shared_ptr<VEC> residual = solver.create_new_vector();
  shared_ptr<VEC> dst = solver.create_new_vector();
  shared_ptr<VEC> src = solver.create_new_vector();

  copy(*src_yy, yy);
  copy(*src_yp, yp);
  copy(*residual, resp);
  copy(*src, b);

  int err = solver.solve_jacobian_system(IDA_mem->ida_tn,
                                         *src_yy,
                                         *src_yp,
                                         *residual,
                                         IDA_mem->ida_cj,
                                         *src,
                                         *dst);
  copy(b, *dst);
  return err;
}


template <typename VEC>
IDAInterface<VEC>::IDAInterface(SundialsInterface<VEC> &bubble) :
  ParameterAcceptor("IDA Solver Parameters"),
  solver(bubble),
  ida_mem(nullptr),
  pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
{
  initial_step_size = 1e-4;
  min_step_size = 1e-6;

  abs_tol = 1e-6;
  rel_tol = 1e-8;
}

template <typename VEC>
IDAInterface<VEC>::~IDAInterface()
{
  if (ida_mem)
    IDAFree(&ida_mem);
}

template <typename VEC>
void IDAInterface<VEC>::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &initial_step_size,
                "Initial step size", "1e-4", Patterns::Double());

  add_parameter(prm, &min_step_size,
                "Min step size", "5e-5", Patterns::Double());

  add_parameter(prm, &abs_tol,
                "Absolute error tolerance", "1e-4", Patterns::Double());

  add_parameter(prm, &rel_tol,
                "Relative error tolerance", "1e-3", Patterns::Double());

  add_parameter(prm, &initial_time,
                "Initial time", "0", Patterns::Double());

  add_parameter(prm, &final_time,
                "Final time", "1", Patterns::Double());

  add_parameter(prm, &outputs_period,
                "Seconds between each output", "1e-1", Patterns::Double());


  add_parameter(prm, &max_order,
                "Maximum order of BDF", "5", Patterns::Integer());


  add_parameter(prm, &max_non_linear_iterations,
                "Maximum number of nonlinear iterations", "10", Patterns::Integer());


  add_parameter(prm, &ignore_algebraic_terms_for_errors,
                "Ignore algebraic terms for error computations", "false",
                Patterns::Bool(),
                "Indicate whether or not to suppress algebraic variables "
                "in the local error test.");

  add_parameter(prm, &ic_type,
                "Initial condition type", "use_y_diff",
                Patterns::Selection("none|use_y_diff|use_y_dot"),
                "This is one of the following thress options for the "
                "initial condition calculation. \n"
                " none: do not try to make initial conditions consistent. \n"
                " use_y_diff: compute the algebraic components of y and differential\n"
                "    components of y_dot, given the differential components of y. \n"
                "    This option requires that the user specifies differential and \n"
                "    algebraic components in the function get_differential_components.\n"
                " use_y_dot: compute all components of y, given y_dot.");

  add_parameter(prm, &reset_type,
                "Initial condition type after restart", "use_y_dot",
                Patterns::Selection("none|use_y_diff|use_y_dot"),
                "This is one of the following thress options for the "
                "initial condition calculation. \n"
                " none: do not try to make initial conditions consistent. \n"
                " use_y_diff: compute the algebraic components of y and differential\n"
                "    components of y_dot, given the differential components of y. \n"
                "    This option requires that the user specifies differential and \n"
                "    algebraic components in the function get_differential_components.\n"
                " use_y_dot: compute all components of y, given y_dot.");

  add_parameter(prm, &ic_alpha,
                "Initial condition Newton parameter", "0.33", Patterns::Double());


  add_parameter(prm, &ic_max_iter,
                "Initial condition Newton max iterations", "5", Patterns::Integer());

  add_parameter(prm, &use_local_tolerances,
                "Use local tolerances", "false", Patterns::Bool());

  add_parameter(prm, &verbose,
                "Show output of time steps", "true", Patterns::Bool());
}


template <typename VEC>
unsigned int IDAInterface<VEC>::start_ode(VEC &solution,
                                          VEC &solution_dot,
                                          const unsigned int max_steps)
{


  AssertThrow(solution.size() == solver.n_dofs(),
              ExcDimensionMismatch(solution.size(), solver.n_dofs()));

  double t = initial_time;
  double h = initial_step_size;
  unsigned int step_number = 0;

  int status;

  // The solution is stored in
  // solution. Here we take only a
  // view of it.

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

  reset_ode(initial_time, solution, solution_dot, initial_step_size, max_steps, true);

  double next_time = 0;

  solver.output_step( 0, solution, solution_dot, 0, initial_step_size);

  while ((t<final_time) && (step_number < max_steps))
    {

      next_time += outputs_period;
      if (verbose)
        {
          pcout << " "//"\r"
                << std::setw(5) << t << " ----> "
                << std::setw(5) << next_time
                << std::endl;
        }
      status = IDASolve(ida_mem, next_time, &t, yy, yp, IDA_NORMAL);

      status = IDAGetLastStep(ida_mem, &h);
      AssertThrow(status == 0, ExcMessage("Error in IDA Solver"));

      copy(solution, yy);
      copy(solution_dot, yp);

      // Check the solution
      bool reset = solver.solver_should_restart(t, step_number, h, solution, solution_dot);


      while (reset)
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

  pcout << std::endl;
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
void IDAInterface<VEC>::reset_ode(double current_time,
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

  if (verbose)
    {
      pcout << "computing consistent initial conditions with the option "
            << type
            << " please be patient."
            << std::endl;
    }

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

  if (verbose)
    {
      pcout << "compute initial conditions: done."
            << std::endl;
    }
}

D2K_NAMESPACE_CLOSE

template class deal2lkit::IDAInterface<BlockVector<double> >;

#ifdef DEAL_II_WITH_MPI

#ifdef DEAL_II_WITH_TRILINOS
template class deal2lkit::IDAInterface<TrilinosWrappers::MPI::Vector>;
template class deal2lkit::IDAInterface<TrilinosWrappers::MPI::BlockVector>;
#endif

#ifdef DEAL_II_WITH_PETSC
template class deal2lkit::IDAInterface<PETScWrappers::MPI::Vector>;
template class deal2lkit::IDAInterface<PETScWrappers::MPI::BlockVector>;
#endif

#endif


#endif
