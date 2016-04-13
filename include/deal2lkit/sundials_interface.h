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

#ifndef _d2k_sundials_interface_h
#define _d2k_sundials_interface_h

#include <deal2lkit/config.h>
#include <deal2lkit/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_view.h>

#ifdef DEAL_II_WITH_MPI
#include "mpi.h"
#endif

#ifdef D2K_WITH_SUNDIALS


using namespace dealii;


D2K_NAMESPACE_OPEN

/**
 * IDA (Implicit Differential-Algebraic solver) is a package contaied in the
 * suite SUNDIALS (SUite of Nonlinear and DIfferential/ALgebraic equation
 * Solvers).  We use this package to provide an interface for problems in the
 * form:
 *
 * @f[ F(t,y,\dot{y})=0, \\ y(t_0)  = y_0 ,   \\ \dot{y}(t_0) = \dot{y}_0.  @f]
 *
 * @remark The integration method used in IDA is the variable-order,
 * variable-coefficient BFD (Backward Differentiation Formula), in
 * fixed-leading-coefficient form.
 *
 * @f[ \sum_{i=0}^{q} \alpha_{n,i}y_{n-1} = h_n \dot{y}_n @f]
 *
 * where @f$ h_n = t_n - t_{n-1}@f$.
 *
 * @remark IDA decompose each vector @f$ y @f$ in a differential and an
 * algebraic part.
 *
 * [see https://computation.llnl.gov/casc/sundials/documentation/ida_guide.pdf
 * for more information about IDA]
 */

template<typename VEC=Vector<double> > class SundialsInterface
{

public :

#ifdef DEAL_II_WITH_MPI
  /**
   * @brief Constructor
   */
  SundialsInterface(const MPI_Comm &communicator=MPI_COMM_WORLD) :
    communicator(communicator) {};

#else
  SundialsInterface() {};
#endif

  /**
   * @brief return the comunicator
   */
  const MPI_Comm &get_comm() const
  {
#ifdef DEAL_II_WITH_MPI
    return communicator;
#else
    return MPI_COMM_WORLD;
#endif
  }


  /**
   * Create a compatible vector with the domain of @f$ F @f$.
   */
  virtual shared_ptr<VEC>
  create_new_vector() const = 0;

  /**
   * Norm used in the solver.
   */
  virtual double
  vector_norm(const VEC &vector) const;

  /**
   * This is a pure virtual function used to get the number of degrees of
   * freedom.
   */
  virtual unsigned int n_dofs() const = 0;

  /**
   * This function is called at the end of each time step of the ODE solver.
   *
   * @remark The conversion between pointers and other forms of vectors need to
   * be done inside the inheriting class.
   *
   * */
  virtual void output_step( const double t,
                            const VEC &solution,
                            const VEC &solution_dot,
                            const unsigned int step_number,
                            const double h) = 0;

  /**
   * This implementation is smart enough to recognize whether the DAE solver
   * should be reset or not.
   *
   * This function returns @c true when the @p solution or @p solution_dot need
   * to be updated (e.g. the number of dofs increases) and @c false otherwise.
   * In few words, this boolean tells to the DEA solver when it has to change
   * its arguments.
   *
   */
  virtual bool solver_should_restart(const double t,
                                     const unsigned int step_number,
                                     const double h,
                                     VEC &solution,
                                     VEC &solution_dot);

  /**
   * Given @p t, @p y, and @p y_dot, this function evaluates the residual of
   * @f[
   *   F(t,y,\dot{y}) = 0
   * @f]
   * and stores the result in the variable @p dst.
   */
  virtual int residual(const double t,
                       const VEC &y,
                       const VEC &y_dot,
                       VEC &dst) = 0;

  /**
   * This function is nothing but the Newton's method used to find the
   * equilibrium point of the (linearized) energy.
   *
   * @p residual is the residual of
   * @p alpha represents the coefficients of the BFD Method
   * @f[
   *   F(t,y,\dot{y}) = 0
   * @f]
   *
   * Given @p t, @p y, and @p y_dot, this function uses the @p residual and
   * @p alpha to solve the linear system:
   * @f[
   *    \mathcal{J}(F) dst = src
   * @f]
   * where @f$ \mathcal{J}(F) @f$ is the last computed Jacobian of @f$ F @f$.
   */
  virtual int solve_jacobian_system(const double t,
                                    const VEC &y,
                                    const VEC &y_dot,
                                    const VEC &residual,
                                    const double alpha,
                                    const VEC &src,
                                    VEC &dst) const = 0;


  /**
   * SundialsInterface uses the Jacobian of @f$ F @f$ to find a minimum for the
   * energy.
   * This function is devoted to the setup of the Jacobian.
   * Setup Jacobian. Compute the Jacobian of the function
   * @f[
   *    \mathcal{J}(F(t, y, \alpha y + c)) = D_y F + alpha D_{\dot{y}} F.
   * @f]
   *
   */
  virtual int setup_jacobian(const double t,
                             const VEC &y,
                             const VEC &y_dot,
                             const VEC &residual,
                             const double alpha) = 0;

  /**
      * compute Jacobian times @p v, and the result is stored in @p dst
      *
      */
  virtual int jacobian_vmult(const VEC &v, VEC &dst) const;

  /**
   * IDA package decomposes variables in differentials and algebraics.
   *
   * This function returns a vector that has 1 on the differential components
   * and 0 on the algebraic ones.
   */
  virtual VEC &differential_components() const;

  /**
   * This function is used to set different tolerances accordingly
   * to the solver.
   */
  virtual VEC &get_local_tolerances() const;

  /**
      * Get the diagonal of the lumped mass matrix and store it in @param diag.
      * @p diag is used as scaling vector for the computation of an approximate
      * L2 norm.
      *
      * By default, @p diag is set to 1.
      */
  virtual void get_lumped_mass_matrix(VEC &diag) const;


#ifdef DEAL_II_WITH_MPI
private:
  /**
   * MPI communicator needed for parallel solver.
   */
  const MPI_Comm &communicator;
#endif


};



D2K_NAMESPACE_CLOSE

#endif

#endif

