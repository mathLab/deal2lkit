//-----------------------------------------------------------
//
//    Copyright (C) 2015-2016 by the deal2lkit authors
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

#ifndef d2k_imex_stepper_h
#define d2k_imex_stepper_h

#include <deal2lkit/config.h>
#include <deal.II/base/config.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>

#ifdef D2K_WITH_SUNDIALS
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/kinsol_interface.h>

#ifdef DEAL_II_WITH_MPI
#include "mpi.h"
#endif

D2K_NAMESPACE_OPEN

/**
 * IMEXStepper solves non-linear time dependent problems with
 * user-defined size of the time step and using Newthon's
 * method for the solution of the non-linear problem.
 * It allows to use the Kinsol solver of \sundials.
 *
 * The user has to provide the following std::functions:
 *  - create_new_vector;
 *  - residual;
 *  - setup_jacobian;
 *  - solve_jacobian_system;
 *  - output_step;
 *  - solver_should_restart;
 *  - get_lumped_mass_matrix (only for kinsol);
 *  - jacobian_vmult (only for kinsol);
 *  - vector_norm (if kinsol is not used).
 *
 */
template<typename VEC=Vector<double> >
class IMEXStepper : public ParameterAcceptor
{
public:

#ifdef DEAL_II_WITH_MPI
  /** Constructor for the IMEXStepper class.
    * Takes a @p name for the section in parameter file
    * and a mpi communicator.
    */
  IMEXStepper(std::string name="",
              MPI_Comm comm = MPI_COMM_WORLD);
#else
  /** Constructor for the IMEXStepper class.
    * Takes a @p name for the section in parameter file.
    */
  IMEXStepper(std::string name="");

#endif

  ~IMEXStepper();

  /** Declare parameters for this class to function properly. */
  virtual void declare_parameters(ParameterHandler &prm);

  /** Evolve. This function returns the final number of steps. */
  unsigned int solve_dae(VEC &solution, VEC &solution_dot);


  /**
   * Compute initial conditions that satisfy \f$ F(y, \dot y, t)=0 \f$.
   */
  void compute_consistent_initial_conditions(const double &t,
                                             VEC &y,
                                             VEC &y_dot);
  /**
   * if initial time is different from final time (i.e.,
   * we are solving a time-dep problem and not a stationay
   * one, return the inverse of dt. If the problem is
   * stationary, returns 0.
   */
  double get_alpha() const;

  /**
   * Set initial time equal to @p t disregarding what
   * is written in the parameter file.
   */
  void set_initial_time(const double &t);

private:

#ifdef DEAL_II_WITH_MPI
  MPI_Comm communicator;
#endif
  /**
   * kinsol solver
   */
  KINSOLInterface<VEC> kinsol;

  void compute_y_dot(const VEC &y, const VEC &prev, const double alpha, VEC &y_dot);

  /** Step size. */
  double step_size;

  /**
   * user defined step_size
   */
  std::string _step_size;

  /**
   * Evaluate step size at time @p t according to the
   * expression stored in _step_size
   */
  double evaluate_step_size(const double &t);

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

  /** Alpha to use in Newton method for the update of the solution. */
  double newton_alpha;

  /** Maximum number of outer iterations for Newton method. */
  unsigned int max_outer_non_linear_iterations;

  /** Maximum number of inner iterations for Newton method. */
  unsigned int max_inner_non_linear_iterations;

  /** Jacobian is updated at each outer iteration and time step */
  bool update_jacobian_continuously;

  /**
   * use kinsol solver true or false
   */
  bool use_kinsol;

  /** Output stream */
  ConditionalOStream pcout;

  /** print useful informations */
  bool verbose;

  /** i max */
  unsigned int n_max_backtracking;

  /** method used for alpha selection*/
  std::string method;

  /**
   *  Line search algorithm with backtracking. The following sequence
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


  /**
   * find solution applying the newton method with given
   * @param t
   * @param alpha
   * @param update_Jacobian
   * @param previous_solution
   * @param solution_dot
   * at the end of the computation, the @p solution_dot is updated as well.
   *
   * this function is called when KINSOL is NOT used
   */
  void  do_newton (const double t,
                   const double alpha,
                   const bool update_Jacobian,
                   const VEC &previous_solution,
                   VEC &solution,
                   VEC &solution_dot);


  /**
   * Compute previous solution from given
   * @param sol
   * @param sol_dot
   * @param alpha
   */
  void compute_previous_solution(const VEC &sol,
                                 const VEC &sol_dot,
                                 const double &alpha,
                                 VEC &prev);

public:

  /**
   * Return a shared_ptr<VEC>. A shared_ptr is needed in order
   * to keep the pointed vector alive, without the need to use a
   * static variable.
   */
  std::function<shared_ptr<VEC>()> create_new_vector;

  /**
   * Compute residual.
   */
  std::function<int(const double t,
                    const VEC &y,
                    const VEC &y_dot,
                    VEC &res)> residual;

  /**
   * Compute Jacobian.
   */
  std::function<int(const double t,
                    const VEC &y,
                    const VEC &y_dot,
                    const double alpha)> setup_jacobian;

  /**
   * Solve linear system.
   */
  std::function<int(const VEC &rhs, VEC &dst)> solve_jacobian_system;

  /**
   * Store solutions to file.
   */
  std::function<void (const double t,
                      const VEC &sol,
                      const VEC &sol_dot,
                      const unsigned int step_number)> output_step;

  /**
   * Evaluate wether the mesh should be refined or not. If so, it
   * refines and interpolate the solutions from the old to the new
   * mesh.
   */
  std::function<bool (const double t,
                      VEC &sol,
                      VEC &sol_dot)> solver_should_restart;

  /**
   * Return the lumped mass matrix vector. It is used
   * by kinsol as scaling factor for the computations of
   * vector's norms.
   */
  std::function<VEC&()> get_lumped_mass_matrix;

  /**
   * Compute the matrix-vector product Jacobian times @p src,
   * and the result is put in @p dst.
   */
  std::function<int(const VEC &src,
                    VEC &dst)> jacobian_vmult;

  /**
   * Return the norm of @p vector. Note that Kinsol uses different
   * norms. By default it returns the l2 norm.
   */
  std::function<double(const VEC &vector)> vector_norm;

private:

  /**
   * Set the std::functions above to trigger an assert if they are not implemented.
   */
  void set_functions_to_trigger_an_assert();

};

D2K_NAMESPACE_CLOSE

#endif


#endif
