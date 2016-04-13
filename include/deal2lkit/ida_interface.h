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
#include <deal.II/base/conditional_ostream.h>

#ifdef D2K_WITH_SUNDIALS


#include <deal2lkit/sundials_interface.h>
#include <deal2lkit/parameter_acceptor.h>



// For time integration.
#include <ida/ida.h>
#include <ida/ida_spils.h>
#include <ida/ida_spgmr.h>
#include <ida/ida_spbcgs.h>
#include <ida/ida_sptfqmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

D2K_NAMESPACE_OPEN

/** Interface to \sundials IDA library.

\dk features an interface for the SUite of Nonlinear and
DIfferential/ALgebraic equation Solvers (\sundials). Each of the
\sundials solvers require the user to derive its class from the base
class SundialsInterface, providing function evaluation, jacobian
evaluation, etc. to the underlying \sundials solver.

The class IDAInterface is a wrapper to the Implicit
Differential-Algebraic solver which is a general purpose solver for
systems of Differential-Algebraic Equations (DAEs).

Citing from the \sundials documentation:

  Consider a system of Differential-Algebraic Equations written in the
  general form

\f[
  \begin{cases}
      F(t,y,\dot y) = 0\, , \\
      y(t_0) = y_0\, , \\
      \dot y (t_0) = \dot y_0\, .
  \end{cases}
\f]

where $y,\dot y$ are vectors in $\R^n$, $t$ is often the time (but can
also be a parametric quantity), and
$F:\R\times\R^n\times\R^n\rightarrow\R^n$. Such problem is solved
using Newton iteration augmented with a line search global
strategy. The integration method used in \ida is the variable-order,
variable-coefficient BDF (Backward Differentiation Formula), in
fixed-leading-coefficient. The method order ranges from 1 to 5, with
the BDF of order $q$ given by the multistep formula

\f[
  \sum\limits_{i=0}^q \alpha_{n,i}\,y_{n-i}=h_n\,\dot y_n\, ,
  \label{eq:bdf}
\f]

where $y_n$ and $\dot y_n$ are the computed approximations of $y(t_n)$ and $\dot y(t_n)$, respectively, and the step size is $h_n=t_n-t_{n-1}$. The coefficients $\alpha_{n,i}$ are uniquely determined by the order $q$, and the history of the step sizes. The application of the BDF method \eqref{eq:bdf} to the DAE system \eqref{eq:dae_system} results in a nonlinear algebraic system to be solved at each time step:

\f[
  G(y_n)\equiv F\left(t_n,y_n,\dfrac{1}{h_n}\sum\limits_{i=0}^q \alpha_{n,i}\,y_{n-i}\right)=0\, .
  \label{eq:nonlinear}
\end{equation}
The Newton method leads to a linear system of the form
\begin{equation}
  J[y_{n(m+1)}-y_{n(m)}]=-G(y_{n(m)})\, ,
  \label{eq:linear}
\f]

where $y_{n(m)}$ is the $m$-th approximation to $y_n$, $J$ is the approximation of the system Jacobian

\f[
  J=\dfrac{\partial G}{\partial y} = \dfrac{\partial F}{\partial y} + \alpha \dfrac{\partial F}{\partial \dot y}\, ,
  \label{eq:jacobian}
\f]

and $\alpha = \alpha_{n,0}/h_n$. It is worthing metioning that the scalar $\alpha$ changes whenever the step size or method order changes.


As far as the solution of the linear system is concerned, the \dk
class SundialsInterface exploits the linear algebra classes included
in the \dealii library (e.g., solvers, preconditioners,
LinearOperator, etc.).

A user that would want to use the \sundials interface, should derive
its problem class from the SundialsInterface class, and
implement all pure virtual methods.
**/
template<typename VEC=Vector<double> >
class IDAInterface : public ParameterAcceptor
{
public:
  /** Constructor for the IDAInterface class. The Solver class is
   * required to have a Solver.solve(VEC &dst, const VEC &src) method
   * that will be called by the time integrator to find out about the
   * solution to a given src. */
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

  /** Output stream */
  ConditionalOStream pcout;
};


D2K_NAMESPACE_CLOSE
#endif


#endif
