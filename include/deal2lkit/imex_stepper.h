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
  unsigned int start_ode(VEC &solution);

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

  /** Maximum number of iterations for Newton method in IC calculation. */
  unsigned int max_non_linear_iterations;
};

D2K_NAMESPACE_CLOSE

#endif


#endif
