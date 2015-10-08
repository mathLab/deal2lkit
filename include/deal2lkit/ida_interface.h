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

#ifndef _d2k_ida_interface_h
#define _d2k_ida_interface_h

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

template<typename VEC=Vector<double> >
class IDAInterface : public ParameterAcceptor
{
public:
  /** Constructor for the IDAInterface class. The Solver class is
   * required to have a Solver.solve(VEC &dst, const
   * VEC &src) method that will be called by the time
   * integrator to find out about the solution to a given src. */
  IDAInterface(SundialsInterface<VEC> &solver);

  /** House cleaning. */
  ~IDAInterface();

  /** Declare parameters for this class to function properly. */
  virtual void declare_parameters(ParameterHandler &prm);

  /** Evolve. This function returns the final number of steps. */
  unsigned int start_ode(VEC &solution,
                         VEC &solution_dot,
                         const unsigned int max_steps);

  /** Clear internal memory, and
  start with clean
  objects. This is useful if
  you need to refine your
  mesh between stesp. */
  void reset_ode(const double t, VEC &y, VEC &yp,
                 double h, unsigned int max_steps,
                 bool first_step);


  /** Final time. */
  double final_time;

  /** Initial time for the ode.*/
  double initial_time;

private:
  /** The bubble membrane poperties. */
  SundialsInterface<VEC> &solver;

  /** Initial step size. */
  double initial_step_size;

  /** Minimum step size. */
  double min_step_size;

  /** Absolute error tolerance for adaptive time stepping. */
  double abs_tol;

  /** Relative error tolerance for adaptive time stepping. */
  double rel_tol;

  /** Maximum number of time steps. */
  unsigned int max_n_steps;

  /** Maximum order of BDF. */
  unsigned int max_order;

  /** Seconds between each output. */
  double outputs_period;

  /** Ignore algebraic terms for errors. */
  bool ignore_algebraic_terms_for_errors;

  /** Type of initial conditions. */
  std::string ic_type;

  /** Type of conditions to be used after a solver restart. */
  std::string reset_type;

  /** Alpha to use in Newton method for IC calculation. */
  double ic_alpha;

  /** Maximum number of iterations for Newton method in IC calculation. */
  unsigned ic_max_iter;

  /** Maximum number of iterations for Newton method during time advancement. */
  unsigned int max_non_linear_iterations;

  /** Initialization flag.*/
  bool is_initialized;

  /** Show the progress of time steps. */
  bool verbose;

  /** Use local tolerances when computing absolute tolerance. */
  bool use_local_tolerances;

  /** Ida memory object. */
  void *ida_mem;

  /** Ida solution vector. */
  N_Vector yy;
  /** Ida solution derivative vector. */
  N_Vector yp;
  /** Ida absolute tolerances vector. */
  N_Vector abs_tolls;
  /** Ida differential components vector. */
  N_Vector diff_id;


};

D2K_NAMESPACE_CLOSE

#endif


#endif
