#ifndef __sak__sundials_interface_h
#define __sak__sundials_interface_h

#include "utilities.h"
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_view.h>

#ifdef DEAL_II_WITH_MPI
#include "mpi.h"

#ifdef DEAL_II_SAK_WITH_SUNDIALS
#include <ida/ida_impl.h>

using namespace dealii;

/**
 * Base class that needs to be inherited by any function that wants
 * to use the DAE time integrator class.
 *
 * Interface for problems of type F(t, y, y_dot) = 0
 *
 * where differential components may be only a part of y.
 */
template<typename VEC=Vector<double> >
class SundialsInterface
{

public :

  SundialsInterface(const MPI_Comm &communicator=MPI_COMM_WORLD) :
    communicator(communicator) {};

  const MPI_Comm &get_comm() const
  {
    return communicator;
  }


  /**
   * Create a compatible vector with the domain of F.
   */
  virtual shared_ptr<VEC>
  create_new_vector() const = 0;

  /** Returns the number of degrees of freedom. Pure virtual function. */
  virtual unsigned int n_dofs() const = 0;

  /** This function is called at the end of each iteration step for
   * the ode solver. Once again, the conversion between pointers and
   * other forms of vectors need to be done inside the inheriting
   * class. */
  virtual void output_step(const double t,
                           const VEC &solution,
                           const VEC &solution_dot,
                           const unsigned int step_number,
                           const double h) = 0;

  /**
   * This function should analyze the solution, and decide wether the dae
   * solver should be restarted.
   *
   * If the function returns true, then both @p solution and @solution_dot may
   * be modified inside this function (also in terms of number of dofs) and the
   * dae solver will use the new ones.
   *
   * Should this function return false, then both @p solution and @p solution_dot
   * are assumed to be unchanged by this function.
   */
  virtual bool solver_should_restart(const double t,
                                     const VEC &solution,
                                     const VEC &solution_dot,
                                     const unsigned int step_number,
                                     const double h);

  /** Compute the residual dst = F(t, y, y_dot). */
  virtual int residual(const double t,
                       const VEC &y,
                       const VEC &sy_dot,
                       VEC &dst) = 0;

  /** Solve the linear system
   *
   * JF dst = src
   *
   * where JF is the last computed Jacobian of F.
   *
   */
  virtual int solve_jacobian_system(const double t,
                                    const VEC &y,
                                    const VEC &y_dot,
                                    const VEC &residual,
                                    const double alpha,
                                    const VEC &src,
                                    VEC &dst) const = 0;


  /** Setup Jacobian. Compute the Jacobian of the function
   *
   * JF = J F(t, y, alpha y + c) = D_y F + alpha D_y_dot F
   *
   */
  virtual int setup_jacobian(const double t,
                             const VEC &y,
                             const VEC &y_dot,
                             const VEC &residual,
                             const double alpha) = 0;

  /**
   * Return a vector with 1 on the differential components, and 0 on the
   * algebraic ones.
   */
  virtual VEC &differential_components() const;

  virtual VEC &get_local_tolerances() const;

private:
  /**
   * MPI communicator for parallel solver.
   */
  const MPI_Comm &communicator;

};

#endif
#endif


#endif
