//-----------------------------------------------------------
//
//    Copyright (C) 2016 by the deal2lkit authors
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

#ifndef _d2k_kinsol_interface_h
#define _d2k_kinsol_interface_h

#include <deal2lkit/config.h>


#ifdef D2K_WITH_SUNDIALS

#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal2lkit/sundials_interface.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/utilities.h>

#ifdef DEAL_II_WITH_MPI
#include "mpi.h"
#endif

#include <kinsol/kinsol_impl.h>

using namespace dealii;

D2K_NAMESPACE_OPEN

/** Interface to \sundials KINSOL library.

\dk features an interface for the SUite of Nonlinear and
DIfferential/ALgebraic equation Solvers (\sundials).

The class KINSOLInterface is a wrapper to the KINSOL solver in \sundials. This is a general-purpose
nonlinear system solver based on Newton-Krylov solver technology.
Developed with KINSOL v.2.8.0.

Citing from the \sundials documentation:

blah blah blah

**/
template<typename VEC=Vector<double> >
class KINSOLInterface : public ParameterAcceptor
{
public:

#ifdef DEAL_II_WITH_MPI

  /**
   * Constructor for the KINSOLInterface class.
   * The first string is the name of the section for
   * ParameterAcceptor
   */
  KINSOLInterface(const std::string name="",
                  const MPI_Comm &mpi_comm = MPI_COMM_WORLD);

  /**
   * Return the comunicator
   */
  const MPI_Comm &get_comm() const ;

#else
  /**
   * Constructor for the KINSOLInterface class.
   * It takes the string containing the name of the section
   * for ParameterHandler
   */
  KINSOLInterface(const std::string name="");

#endif

  /**
   * House cleaning.
   */
  ~KINSOLInterface();

  /**
   * Declare parameters for this class to function properly.
   */
  virtual void declare_parameters( ParameterHandler &prm );

  /**
   * Initializes the solver with the initial guess and the residual function.
   * The system size is set inside this function. So, when a mesh refinement is done,
   * you need to call this function in order to make kinsol aware of the changes.
   */
  void initialize_solver( VEC &initial_guess );

  /**
   * Set the scaling for the solver.
   * @p uscale is the scaling vector for the solution and @p fscale is the scaling for the Jacobian.
   * If this function is not called by the user, a scaling factor equal to 1 is assumed.
         */
  void set_scaling_vectors( const VEC &uscale, const VEC &fscale );

  /**
   * Set the constraint.
   *  Citing from the \sundials documentation:
   *  if constraint[i] is
   *     0.0 then no constraint is imposed on u_i.
   *     1.0 then u_i will be constrained to be u_i >= 0.0.
   *    −1.0 then u_i will be constrained to be u_i <= 0.0.
   *     2.0 then u_i will be constrained to be u_i > 0.0.
   *    −2.0 then u_i will be constrained to be u_i < 0.0.
   */
  void set_constraint_vector( const VEC &constraint );

  /**
   * solve the non-linear system
        */
  int solve( VEC &solution );

  /**
   * this function has to be implemented by the user
   * and it must return a shared pointer to a VEC vector.
   */
  std::function<shared_ptr<VEC>()> create_new_vector;

  /** standard function computing residuals */
  std::function<int(const VEC &y, VEC &res)> residual;

  /** standard function computing the Jacobian */
  std::function<int(const VEC &y)> setup_jacobian;

  /** standard function solving linear system */
  std::function<int(const VEC &res, VEC &dst)> solve_linear_system;

  /** standard function multiplying the Jacobian to a vector */
  std::function<int(const VEC &v, VEC &dst )> jacobian_vmult;

private:

  /**
   * This function is executed at construction time to set the std::function above to trigger an assert if they are not implemented.
   */
  void set_functions_to_trigger_an_assert();

  /** Strategy used by the solver:
    *    newton        = basic Newton iteration
    *    global_newton = Newton with line search
    *    fixed_point   = fixed-point iteration with Anderson Acceleration
    *    picard        = Picard iteration with Anderson Acceleration
    */
  std::string strategy;

  /** Maximum number of iterations */
  unsigned int max_iterations;

  /**
   * Input scalar tolerance for residuals.
   *  This define the condition (small residual) for a successful completion of KINSOL
   */
  double ftol;

  /**
   * The Newton method will terminate when the scaled norm of the update is smaller than steptol
   */
  double steptol;

  /**
   * Maximum number of nonlinear iteration that can be done without a Jacobian update
   */
  double mbset;

  /**
   * dimension of the system
   */
  unsigned int system_size;

  /**
   * Initialization flag.
   */
  bool is_initialized;

  /**
   * scaling flag: if the scaling have been provided by the user
   * calling the function set_scaling_vectors();
   */
  bool scaling_is_set;

  /**
   * Level of verbosity of kinsol solver 0,1,2,3
   */
  unsigned int verbosity;

  /**
   * if true the internal direct solver of KINSOL is used.
   * note that the direct solver of KINSOL work only in serial.
   * if false the solver provided by the user is used.
   */
  bool use_internal_solver;

  /**
   * KINSOL solution vector.
   * After initialize_solver is called this contains the initial guess.
   * After solve is called this contains the solution to the equations.
   */
  N_Vector solution;

  /**
   * Internal scaling vector for the solution
   */
  N_Vector u_scale;

  /**
   * Internal scaling vector for the Jacobian
   */
  N_Vector f_scale;

#ifdef DEAL_II_WITH_MPI
  /**
   * local size of the system
   */
  unsigned int local_system_size;

  /**
   * MPI communicator needed for parallel solver.
   */
  const MPI_Comm &communicator;
#endif

  /**
   * Output stream
   */
  ConditionalOStream pcout;

  /**
   * Kinsol memory object.
   */
  void *kin_mem;


};

D2K_NAMESPACE_CLOSE

#endif

#endif
