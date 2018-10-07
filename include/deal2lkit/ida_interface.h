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
//---------------------------------------------------------------

#ifndef d2k_ida_interface_h
#define d2k_ida_interface_h

#include <deal.II/base/config.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_view.h>

#include <deal2lkit/config.h>
#include <deal2lkit/utilities.h>

#ifdef D2K_WITH_SUNDIALS


#  include <deal2lkit/parameter_acceptor.h>
#  include <ida/ida.h>
#  include <ida/ida_spbcgs.h>
#  include <ida/ida_spgmr.h>
#  include <ida/ida_spils.h>
#  include <ida/ida_sptfqmr.h>
#  include <nvector/nvector_serial.h>
#  include <sundials/sundials_math.h>
#  include <sundials/sundials_types.h>

#  ifdef DEAL_II_WITH_MPI
#    include "mpi.h"
#  endif

D2K_NAMESPACE_OPEN

/** Interface to \sundials IDA library.
 *
 * \dk features an interface for the SUite of Nonlinear and
 * DIfferential/ALgebraic equation Solvers (\sundials).
 *
 * The class IDAInterface is a wrapper to the Implicit
 * Differential-Algebraic solver which is a general purpose solver for
 * systems of Differential-Algebraic Equations (DAEs).
 *
 * The user has to provide the implmentation of the following std::functions:
 *  - create_new_vector;
 *  - residual;
 *  - setup_jacobian;
 *  - solve_jacobian_system;
 *  - output_step;
 *  - solver_should_restart;
 *  - differential_components.
 *
 * Citing from the \sundials documentation:
 *
 *   Consider a system of Differential-Algebraic Equations written in the
 *   general form
 *
 * \f[
 *   \begin{cases}
 *       F(t,y,\dot y) = 0\, , \\
 *       y(t_0) = y_0\, , \\
 *       \dot y (t_0) = \dot y_0\, .
 *   \end{cases}
 * \f]
 *
 * where \f$y,\dot y\f$ are vectors in \f$\R^n\f$, \f$t\f$ is often the time
 * (but can also be a parametric quantity), and
 * \f$F:\R\times\R^n\times\R^n\rightarrow\R^n\f$. Such problem is solved
 * using Newton iteration augmented with a line search global
 * strategy. The integration method used in ida is the variable-order,
 * variable-coefficient BDF (Backward Differentiation Formula), in
 * fixed-leading-coefficient. The method order ranges from 1 to 5, with
 * the BDF of order \f$q\f$ given by the multistep formula
 *
 * \f[
 *   \sum\limits_{i=0}^q \alpha_{n,i}\,y_{n-i}=h_n\,\dot y_n\, ,
 *   \label{eq:bdf}
 * \f]
 *
 * where \f$y_n\f$ and \f$\dot y_n\f$ are the computed approximations of
 * \f$y(t_n)\f$ and \f$\dot y(t_n)\f$, respectively, and the step size is
 * \f$h_n=t_n-t_{n-1}\f$. The coefficients \f$\alpha_{n,i}\f$ are uniquely
 * determined by the order \f$q\f$, and the history of the step sizes. The
 * application of the BDF method to the DAE system results in a nonlinear
 * algebraic system to be solved at each time step:
 *
 * \f[
 *   G(y_n)\equiv F\left(t_n,y_n,\dfrac{1}{h_n}\sum\limits_{i=0}^q
 * \alpha_{n,i}\,y_{n-i}\right)=0\, . \label{eq:nonlinear} \end{equation} The
 * Newton method leads to a linear system of the form \begin{equation}
 *   J[y_{n(m+1)}-y_{n(m)}]=-G(y_{n(m)})\, ,
 *   \label{eq:linear}
 * \f]
 *
 * where \f$y_{n(m)}\f$ is the \f$m\f$-th approximation to \f$y_n\f$, \f$J\f$ is
 * the approximation of the system Jacobian
 *
 * \f[
 *   J=\dfrac{\partial G}{\partial y} = \dfrac{\partial F}{\partial y} + \alpha
 * \dfrac{\partial F}{\partial \dot y}\, , \label{eq:jacobian} \f]
 *
 * and \f$\alpha = \alpha_{n,0}/h_n\f$. It is worthing metioning that the
 * scalar \f$\alpha\f$ changes whenever the step size or method order
 * changes.
 *
 *
 * As far as the solution of the linear system is concerned, the \dk
 * class SundialsInterface exploits the linear algebra classes included
 * in the \dealii library (e.g., solvers, preconditioners,
 * LinearOperator, etc.).
 *
 * A user that would want to use the \sundials interface, should derive
 * its problem class from the SundialsInterface class, and
 * implement all pure virtual methods.
 */
template <typename VEC = Vector<double>>
class IDAInterface : public ParameterAcceptor
{
public:
#  ifdef DEAL_II_WITH_MPI
  /**
   * Constructor for the IDAInterface class.
   */
  IDAInterface(const std::string name     = "",
               const MPI_Comm    mpi_comm = MPI_COMM_WORLD);
#  else
  /**
   * Constructor for the IDAInterface class.
   */
  IDAInterface(const std::string name = "");
#  endif

  /** House cleaning. */
  ~IDAInterface();

  /** Declare parameters for this class to function properly. */
  virtual void declare_parameters(ParameterHandler &prm);

  /** Evolve. This function returns the final number of steps. */
  unsigned int solve_dae(VEC &solution, VEC &solution_dot);

  /**
   * Clear internal memory, and start with clean objects. This
   * function is called when the simulation start and when the mesh is
   * refined.
   */
  void reset_dae(const double t, VEC &y, VEC &yp, double h, bool first_step);

  /**
   * Return a shared_ptr<VEC>. A shared_ptr is needed in order
   * to keep the pointed vector alive, without the need to use a
   * static variable.
   */
  std::function<shared_ptr<VEC>()> create_new_vector;

  /**
   * Compute residual.
   */
  std::function<int(const double t, const VEC &y, const VEC &y_dot, VEC &res)>
    residual;

  /**
   * Compute Jacobian.
   */
  std::function<
    int(const double t, const VEC &y, const VEC &y_dot, const double alpha)>
    setup_jacobian;

  /**
   * Solve linear system.
   */
  std::function<int(const VEC &rhs, VEC &dst)> solve_jacobian_system;

  /**
   * Store solutions to file.
   */
  std::function<void(const double       t,
                     const VEC &        sol,
                     const VEC &        sol_dot,
                     const unsigned int step_number)>
    output_step;

  /**
   * Evaluate wether the mesh should be refined or not. If so, it
   * refines and interpolate the solutions from the old to the new
   * mesh.
   */
  std::function<bool(const double t, VEC &sol, VEC &sol_dot)>
    solver_should_restart;

  /**
   * Return a vector whose component are 1 if the corresponding
   * dof is differential, 0 if algebraic.
   */
  std::function<VEC &()> differential_components;

  /**
   * Return a vector whose components are the weights used by IDA to
   * compute the vector norm. The implementation of this function
   * is optional.
   */
  std::function<VEC &()> get_local_tolerances;



  /**
   * Set initial time equal to @p t disregarding what is written
   * in the parameter file.
   */
  void set_initial_time(const double &t);

private:
  /**
   * This function is executed at construction time to set the
   * std::function above to trigger an assert if they are not
   * implemented.
   */
  void set_functions_to_trigger_an_assert();


  /** Final time. */
  double final_time;

  /** Initial time for the ode.*/
  double initial_time;

  /** Initial step size. */
  double initial_step_size;

  /** Minimum step size. */
  double min_step_size;

  /** Absolute error tolerance for adaptive time stepping. */
  double abs_tol;

  /** Relative error tolerance for adaptive time stepping. */
  double rel_tol;

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

#  ifdef DEAL_II_WITH_MPI
  MPI_Comm communicator;
#  endif

  /** Output stream */
  ConditionalOStream pcout;

  unsigned int system_size;

  unsigned int local_system_size;
};


D2K_NAMESPACE_CLOSE
#endif


#endif
