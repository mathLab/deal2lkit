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

#ifndef _d2k_imex_stepper_h
#define _d2k_imex_stepper_h

#include <deal2lkit/config.h>
#include <deal.II/base/config.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>

#ifdef D2K_WITH_SUNDIALS
// For time integration.
#include <ida/ida.h>
#include <ida/ida_spils.h>
#include <ida/ida_spgmr.h>
#include <ida/ida_spbcgs.h>
#include <ida/ida_sptfqmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>


#include <deal2lkit/sundials_interface.h>
#include <deal2lkit/parameter_acceptor.h>


D2K_NAMESPACE_OPEN

/**
 *
 */
template<typename VEC=Vector<double> >
class IMEXStepper : public ParameterAcceptor
{
public:
  /** Constructor for the IDAInterface class. The Solver class is
   * required to have a Solver.solve(VEC &dst, const VEC &src) method
   * that will be called by the time integrator to find out about the
   * solution to a given src. */
  IMEXStepper(SundialsInterface<VEC> &solver,
              const double &step_size=1e-3,
              const double &initial_time=0.0,
              const double &final_time=1.0
             );

  /** Declare parameters for this class to function properly. */
  virtual void declare_parameters(ParameterHandler &prm);

  /** Evolve. This function returns the final number of steps. */
  unsigned int start_ode(VEC &solution, VEC &solution_dot);

  /** Clear internal memory, and
   start with clean
   objects. This is useful if
   you need to refine your
   mesh between stesp. */
  void reset_ode(const double t, VEC &y,
                 double h, unsigned int max_steps,
                 bool first_step);

private:
  /** The solver interface. */
  SundialsInterface<VEC> &interface;

  /** Step size. */
  double step_size;

  /** Initial time for the ode.*/
  double initial_time;

  /** Final time. */
  double final_time;

  /** Absolute error tolerance for non linear iterations. */
  double abs_tol;

  /** Relative error tolerance for non linear iterations. */
  double rel_tol;

  /** Seconds between each output. */
  unsigned int output_period;

  /** Alpha to use in Newton method for IC calculation. */
  double newton_alpha;

  /** Maximum number of outer iterations for Newton method. */
  unsigned int max_outer_non_linear_iterations;

  /** Maximum number of inner iterations for Newton method. */
  unsigned int max_inner_non_linear_iterations;

  /** Jacobian is updated at each outer iteration and time step */
  bool update_jacobian_continuously;

  /** Output stream */
  ConditionalOStream pcout;

  /** print useful informations */
  bool verbose;

  /** i max */
  unsigned int n_max_backtracking;

  /** method used for alpha selection*/
  std::string method;

  /**
   *  line search algorithm with backtracking.The following sequence
   *  of Newton relaxation parameters is tested: 1, 1/2, 1/4,...,2^-i.
   *  it returns the selected alpha and solution, solution_dot and
   *  residual are accordingly updated.
   */
  double line_search_with_backtracking(const VEC &update,
                                       const VEC &prev_sol,
                                       const double &alpha,
                                       const double &t,
                                       VEC &sol,
                                       VEC &sol_dot,
                                       VEC &residual);

};

D2K_NAMESPACE_CLOSE

#endif


#endif
