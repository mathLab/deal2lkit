#ifdef DEAL_II_SAK_WITH_SUNDIALS

#include "dae_time_integrator.h"
#include "sundials_interface.h"

#include <deal.II/base/utilities.h>
#include <deal.II/lac/block_vector.h>
#ifdef DEAL_II_WITH_TRILINOS
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#endif

#include <iostream>

#include <nvector/nvector_parallel.h>
#include <nvector/nvector_serial.h>

using namespace dealii;
using namespace std;

void copy(TrilinosWrappers::MPI::Vector &dst, const N_Vector &src)
{
  IndexSet is = dst.locally_owned_elements();
  AssertDimension(is.n_elements(), NV_LOCLENGTH_P(src));
  for (unsigned int i=0; i<is.n_elements(); ++i)
    {
      dst[is.nth_index_in_set(i)] = NV_Ith_P(src, i);
    }
}

void copy(N_Vector &dst, const TrilinosWrappers::MPI::Vector &src)
{
  IndexSet is = src.locally_owned_elements();
  AssertDimension(is.n_elements(), NV_LOCLENGTH_P(dst));
  for (unsigned int i=0; i<is.n_elements(); ++i)
    {
      NV_Ith_P(dst, i) = src[is.nth_index_in_set(i)];
    }
}

void copy(TrilinosWrappers::MPI::BlockVector &dst, const N_Vector &src)
{
  IndexSet is = dst.locally_owned_elements();
  AssertDimension(is.n_elements(), NV_LOCLENGTH_P(src));
  for (unsigned int i=0; i<is.n_elements(); ++i)
    {
      dst[is.nth_index_in_set(i)] = NV_Ith_P(src, i);
    }
}

void copy(N_Vector &dst, const TrilinosWrappers::MPI::BlockVector &src)
{
  IndexSet is = src.locally_owned_elements();
  AssertDimension(is.n_elements(), NV_LOCLENGTH_P(dst));
  for (unsigned int i=0; i<is.n_elements(); ++i)
    {
      NV_Ith_P(dst, i) = src[is.nth_index_in_set(i)];
    }
}


void copy(BlockVector<double> &dst, const N_Vector &src)
{
  AssertDimension((unsigned int)NV_LOCLENGTH_P(src), dst.size());
  for (unsigned int i=0; i<dst.size(); ++i)
    {
      dst[i] = NV_Ith_P(src, i);
    }
}

void copy(N_Vector &dst, const BlockVector<double> &src)
{
  AssertDimension((unsigned int)NV_LOCLENGTH_P(dst), src.size());
  for (unsigned int i=0; i<src.size(); ++i)
    {
      NV_Ith_P(dst, i) = src[i];
    }
}




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
DAETimeIntegrator<VEC>::DAETimeIntegrator(SundialsInterface<VEC> &bubble) :
  ParameterAcceptor("IDA Solver Parameters"),
  solver(bubble),
  is_initialized(false)
{
  initial_step_size = 1e-4;
  min_step_size = 1e-6;

  abs_tol = 1e-6;
  rel_tol = 1e-8;

  ida_mem = IDACreate();
  is_initialized = true;

}

template <typename VEC>
DAETimeIntegrator<VEC>::~DAETimeIntegrator()
{
  if (ida_mem)
    IDAFree(&ida_mem);
}

template <typename VEC>
void DAETimeIntegrator<VEC>::declare_parameters(ParameterHandler &prm)
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

  add_parameter(prm, &ic_alpha,
                "Initial condition Newton parameter", "0.33", Patterns::Double());


  add_parameter(prm, &ic_max_iter,
                "Initial condition Newton max iterations", "5", Patterns::Integer());

  add_parameter(prm, &use_local_tolerances,
                "Use local tolerances", "false", Patterns::Bool());
}


template <typename VEC>
unsigned int DAETimeIntegrator<VEC>::start_ode(VEC &solution,
                                               VEC &solution_dot,
                                               const unsigned int max_steps)
{


  AssertThrow(solution.size() == solver.n_dofs(),
              ExcDimensionMismatch(solution.size(), solver.n_dofs()));

  AssertThrow(is_initialized, ExcMessage("Not Initialized!"));

  double t = initial_time;
  double h = initial_step_size;
  unsigned int step_number = 0;

  int status;

  // The solution is stored in
  // solution. Here we take only a
  // view of it.

  IndexSet is = solution.locally_owned_elements();


  yy        = N_VNew_Parallel(solver.get_comm(), is.n_elements(), solver.n_dofs());
  yp        = N_VNew_Parallel(solver.get_comm(), is.n_elements(), solver.n_dofs());
  diff_id   = N_VNew_Parallel(solver.get_comm(), is.n_elements(), solver.n_dofs());
  abs_tolls = N_VNew_Parallel(solver.get_comm(), is.n_elements(), solver.n_dofs());

  reset_ode(initial_time, solution, solution_dot, initial_step_size, max_steps);

  double next_time = 0;

  solver.output_step( 0, solution, solution_dot, 0, initial_step_size);

  while ((t<final_time) && (step_number < max_steps))
    {

      next_time += outputs_period;
      cout << "\r"
           << setw(5) << t << " ----> "
           << setw(5) << next_time
           << " ### ";
      status = IDASolve(ida_mem, next_time, &t, yy, yp, IDA_NORMAL);

      status = IDAGetLastStep(ida_mem, &h);
      AssertThrow(status == 0, ExcMessage("Error in IDA Solver"));
      cout << setw(4) << "Step " << step_number
           << setw(4) << ", t = " << t
           << setw(4) << ", h = " << h;

      copy(solution, yy);
      copy(solution_dot, yp);

      // Check the solution
      bool reset = solver.solver_should_restart(t, solution, solution_dot, step_number, h);


      solver.output_step(t, solution, solution_dot,  step_number, h);

      if ( reset == true )
        {
          double frac = 0;
          int k = 0;
          IDAGetLastOrder(ida_mem, &k);
          frac = std::pow((double)k,2.);
          reset_ode(t, solution, solution_dot,
                    h/frac, max_steps);
        }


      step_number++;
    }

  cout << endl;
  // Free the vectors which are no longer used.
  N_VDestroy_Parallel(yy);
  N_VDestroy_Parallel(yp);
  N_VDestroy_Parallel(abs_tolls);
  N_VDestroy_Parallel(diff_id);

  return step_number;
}

template <typename VEC>
void DAETimeIntegrator<VEC>::reset_ode(double current_time,
                                       VEC &solution,
                                       VEC &solution_dot,
                                       double current_time_step,
                                       unsigned int max_steps)
{
  if (ida_mem)
    IDAFree(&ida_mem);

  ida_mem = IDACreate();


  // Free the vectors which are no longer used.
  if (yy)
    {
      N_VDestroy_Parallel(yy);
      N_VDestroy_Parallel(yp);
      N_VDestroy_Parallel(abs_tolls);
      N_VDestroy_Parallel(diff_id);
    }

  int status;
  Assert(solution.size() == solver.n_dofs(),
         ExcDimensionMismatch(solution.size(), solver.n_dofs()));

  Assert(solution_dot.size() == solver.n_dofs(),
         ExcDimensionMismatch(solution_dot.size(), solver.n_dofs()));


  IndexSet is = solution.locally_owned_elements();

  yy        = N_VNew_Parallel(solver.get_comm(), is.n_elements(), solver.n_dofs());
  yp        = N_VNew_Parallel(solver.get_comm(), is.n_elements(), solver.n_dofs());
  diff_id   = N_VNew_Parallel(solver.get_comm(), is.n_elements(), solver.n_dofs());
  abs_tolls = N_VNew_Parallel(solver.get_comm(), is.n_elements(), solver.n_dofs());

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

  if (ic_type == "use_y_dot")
    {
      // (re)initialization of the vectors
      IDACalcIC(ida_mem, IDA_Y_INIT, current_time+current_time_step);
      IDAGetConsistentIC(ida_mem, yy, yp);

      copy(solution, yy);
      copy(solution_dot, yp);
    }
  else if (ic_type == "use_y_diff")
    {
      IDACalcIC(ida_mem, IDA_YA_YDP_INIT, current_time+current_time_step);
      IDAGetConsistentIC(ida_mem, yy, yp);

      copy(solution, yy);
      copy(solution_dot, yp);
    }
}


#ifdef DEAL_II_WITH_TRILINOS
template class DAETimeIntegrator<TrilinosWrappers::MPI::Vector>;
template class DAETimeIntegrator<TrilinosWrappers::MPI::BlockVector>;
#endif

template class DAETimeIntegrator<BlockVector<double> >;
#endif
