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

#ifndef _d2k_arkode_interface_h
#define _d2k_arkode_interface_h

#include <deal2lkit/config.h>


#ifdef D2K_WITH_SUNDIALS

#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/utilities.h>

#ifdef DEAL_II_WITH_MPI
#include "mpi.h"
#endif

#include <arkode/arkode_impl.h>

using namespace dealii;

D2K_NAMESPACE_OPEN

template<typename VEC=Vector<double> >
class ARKodeInterface : public ParameterAcceptor
{
public:

#ifdef DEAL_II_WITH_MPI
  ARKodeInterface(const std::string name="",
                  const MPI_Comm = MPI_COMM_WORLD);
#else
  ARKodeInterface(const std::string name="");
#endif

  ~ARKodeInterface();

  /**
   * Initialize the solver: set @p t0 as initial time
   * and @p y0 as the initial condition vector.
   * Internally, it sets the functions defining the
   * explicit and/or implici portion of the right hand side.
   */
  void initialize(const double &t0,
                  const VEC &y0);

  void solve(VEC &solution);

  std::function<shared_ptr<VEC>() > create_new_vector;

  std::function<int(const double &t,
                    const VEC &y,
                    VEC &expl_rhs)> explicit_rhs;

  std::function<int(const double &t,
                    const VEC &y,
                    VEC &impl_rhs)> implicit_rhs;

  std::function<int(const double &t)> mass_matrix;

  std::function<int(const VEC &src,
                    VEC &dst)> mass_matrix_vmult;


  std::function<int(const double &gamma,
                    const VEC &src,
                    VEC &dst)> solve_linear_system;

  std::function<int(const VEC &src,
                    VEC &dst)> solve_mass_system;


  std::function<int(const double &t, const VEC &y)> setup_jacobian;

  /**
   * This function is called by ARKode to resize
   * internal vectors.
   */
  std::function<int(const VEC &exemplar,
                    VEC &v)> resize_vector;





  /**
   * Declare parameters for this class to function properly.
   */
  virtual void declare_parameters( ParameterHandler &prm );

private:

#ifdef DEAL_II_WITH_MPI
  MPI_Comm communicator;
  /**
   * local size of the system
   */
  unsigned int local_system_size;
#endif
  unsigned int system_size;


  /**
   * This function is executed at construction time to set the
   * std::function above to trigger an assert if they are not
   * implemented.
   */
  void set_functions_to_trigger_an_assert();

  N_Vector internal_solution;

  /**
   * Output stream
   */
  ConditionalOStream pcout;

  /**
   * ARKode memory object.
   */
  void *ark_mem;

  /**
   * relative tolerance
   */
  double rtol;

  /**
   * absolute tolerance
   */
  double atol;

  /**
   * initial time
   */
  double itime;

  /**
   * final time
   */
  double ftime;

};


D2K_NAMESPACE_CLOSE
#endif //D2K_WITH_SUNDIALS
#endif
