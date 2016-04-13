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

#include <deal2lkit/kinsol_interface.h>
#ifdef D2K_WITH_SUNDIALS
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/generic_linear_algebra.h>

#ifdef DEAL_II_WITH_PETSC
#include <deal.II/lac/petsc_solver.h>
#ifdef DEAL_II_WITH_MPI
#include <deal.II/lac/petsc_parallel_vector.h>
#endif
#endif

#ifdef DEAL_II_WITH_TRILINOS
#include <deal.II/lac/trilinos_block_vector.h>
#ifdef DEAL_II_WITH_MPI
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#endif
#endif

#include <deal2lkit/utilities.h>
#include <kinsol/kinsol_dense.h>

D2K_NAMESPACE_OPEN

// protect the helper functions
namespace
{

  /** helper function to interface residual to sundials */
  template<typename VEC>
  int kinsol_residual( N_Vector y, N_Vector res, void *user_data)
  {

    KINSOLInterface<VEC> &solver = *static_cast<KINSOLInterface<VEC> *>(user_data);


    shared_ptr<VEC> loc_y = solver.create_new_vector();
    shared_ptr<VEC> loc_res = solver.create_new_vector();

    copy(*loc_y, y);
    copy(*loc_res, res);

    int ret = solver.residual(*loc_y, *loc_res);

    copy(res, *loc_res);
    return ret;
  }

  /** helper function to interface setup_jacobian to sundials */
  template<typename VEC>
  int kinsol_setup_jacobian( KINMem kin_mem )
  {
    KINSOLInterface<VEC> &solver = *static_cast<KINSOLInterface<VEC> *> (kin_mem->kin_user_data);

    shared_ptr<VEC> loc_y = solver.create_new_vector();

    copy( *loc_y, kin_mem->kin_uu );
    int err = solver.setup_jacobian( *loc_y );
    return err;
  }

  /** helper function to interface solve_linear_system to sundials */
  template<typename VEC>
  int kinsol_solve_linear_system(KINMem kin_mem,
                                 N_Vector x,
                                 N_Vector b,
                                 double   *sJpnorm,
                                 double   *sFdotJp )
  {
    KINSOLInterface<VEC> &solver = *static_cast<KINSOLInterface<VEC> *> (kin_mem->kin_user_data);

    //shared_ptr<VEC> loc_y   = solver.create_new_vector();
    shared_ptr<VEC> loc_res = solver.create_new_vector();
    shared_ptr<VEC> dst     = solver.create_new_vector();
    shared_ptr<VEC> loc_b   = solver.create_new_vector();

    //copy( *loc_y, x );
    copy( *loc_res, kin_mem->kin_fval);

    int err = solver.solve_linear_system(*loc_res, *dst );

    copy(x, *dst);
    err += solver.jacobian_vmult( *dst, *loc_b );

    copy(b, *loc_b);

    *sJpnorm = N_VWL2Norm(b, kin_mem->kin_fscale);
    N_VProd(b, kin_mem->kin_fscale, b);
    N_VProd(b, kin_mem->kin_fscale, b);
    *sFdotJp = N_VDotProd(kin_mem->kin_fval, b);
    return err;

  }

}

template <typename VEC>
void KINSOLInterface<VEC>::set_functions_to_trigger_an_assert()
{

  create_new_vector = []() ->shared_ptr<VEC>
  {
    shared_ptr<VEC> p;
    AssertThrow(false, ExcPureFunctionCalled());
    return p;
  };

  residual = [](const VEC &, VEC &) ->int
  {
    int ret=0;
    AssertThrow(false, ExcPureFunctionCalled());
    return ret;
  };

  setup_jacobian = [](const VEC &) ->int
  {
    int ret=0;
    AssertThrow(false, ExcPureFunctionCalled());
    return ret;
  };

  solve_linear_system = [](const VEC &, VEC &) ->int
  {
    int ret=0;
    AssertThrow(false, ExcPureFunctionCalled());
    return ret;
  };

  jacobian_vmult = [](const VEC &, VEC &) ->int
  {
    int ret =0;
    AssertThrow(false, ExcPureFunctionCalled());
    return ret;
  };

}

// Class constructor:
#ifdef DEAL_II_WITH_MPI


template <typename VEC>
KINSOLInterface<VEC>::KINSOLInterface(const std::string name,
                                      const MPI_Comm &mpi_comm):
  ParameterAcceptor(name),
  is_initialized(false),
  scaling_is_set(false),
  communicator(mpi_comm),
  pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_comm)==0),
  kin_mem(nullptr)
{
  set_functions_to_trigger_an_assert();
}

#else

template <typename VEC>
KINSOLInterface<VEC>::KINSOLInterface(const std::string name):
  ParameterAcceptor(name),
  is_initialized(false),
  pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
  kin_mem(nullptr)
{
  set_functions_to_trigger_an_assert();
}

#endif

// Class destructor:
template <typename VEC>
KINSOLInterface<VEC>::~KINSOLInterface()
{
  if (kin_mem)
    KINFree(&kin_mem);
}

// Parameters parsing and initialization of the parameters:
template <typename VEC>
void KINSOLInterface<VEC>::declare_parameters(ParameterHandler &prm)
{

  add_parameter(prm, &max_iterations,
                "Maximum number of iterations","200",
                Patterns::Integer(0),
                "The non-linear solver will stop when reaching this number of"
                "Newton iterations no matter what.");

  add_parameter(prm, &ftol,
                "Tolerance for residuals", "1e-9",
                Patterns::Double(0),
                "This define the condition (small residual) for a successful completion of KINSOL.\n"
                " 0 means the default value for KINSOL");

  add_parameter(prm, &steptol,
                "Step tolerance", "1e-11",
                Patterns::Double(0),
                "The Newton method will terminate when the maximum scaled step is below the given tolerance." );

  add_parameter(prm, &verbosity,
                "Level of verbosity of the KINSOL solver", "0", Patterns::Integer(0,3),
                "Allowed values [0,3]");

  add_parameter(prm, &strategy,"Strategy" ,"newton",
                Patterns::Selection("newton|global_newton|fixed_point|picard"),
                "newton        = basic Newton iteration \n"
                "global_newton = Newton with line search \n"
                "fixed_point   = fixed-point iteration with Anderson Acceleration \n"
                "picard        = Picard iteration with Anderson Acceleration");

  add_parameter(prm, &mbset,
                "Maximum number of iteration before Jacobian update", "10",
                Patterns::Integer(1),
                "Maximum number of nonlinear iterations that can be done with an outdated Jacobian.\n"
                "If set to 1 the Jacobian is updated at each nonlinear iteration");

  add_parameter(prm, &use_internal_solver,
                "Use internal KINSOL direct solver", "false",
                Patterns::Bool(),
                "If true the dense direct linear solver of KINSOL is used");

}

// Initialization of the solver:
template <typename VEC>
void KINSOLInterface<VEC>::initialize_solver( VEC &initial_guess )
{
  int status;

  // check if the solver is initialized. If it is reset:
  if (kin_mem)
    {
      KINFree(&kin_mem);
      is_initialized = false;
      scaling_is_set = false;
    }
  // create the KINSOL memory object:
  kin_mem = KINCreate();

  // get the size of the initial guess as the size of the system:
  system_size    = initial_guess.size();
  // need to feed KINSOL init with a N_VECTOR. For simplicity the initial guess is stored in the solution vector.

#ifdef DEAL_II_WITH_MPI
  IndexSet is = initial_guess.locally_owned_elements();
  local_system_size = is.n_elements();
  solution    = N_VNew_Parallel(communicator, local_system_size, system_size);
#else
  solution    = N_VNew_Serial(system_size);
#endif

  // copy the initial guess on the solution:
  copy( solution, initial_guess );

  // pass the pointer to the class to kinsol:
  status = KINSetUserData(kin_mem, (void *) this);
  Assert(status == KIN_SUCCESS, ExcMessage("Error initializing KINSOL. KINInit failed."));

  //   1- maximum iterations:
  status = KINSetNumMaxIters(kin_mem, max_iterations );
  Assert(status == KIN_SUCCESS, ExcMessage("Error initializing KINSOL. KINSetNumMaxIters failed."));

  //   2- ftol:
  status = KINSetFuncNormTol(kin_mem, ftol);
  Assert(status == KIN_SUCCESS, ExcMessage("Error initializing KINSOL. KINSetFuncNormTol failed."));

  //   3- steptol:
  status = KINSetScaledStepTol(kin_mem, steptol);
  Assert(status == KIN_SUCCESS, ExcMessage("Error initializing KINSOL. KINSetFuncNormTol failed."));

  // initialize with the helper function:
  status = KINInit(kin_mem, kinsol_residual<VEC> , solution);
  Assert(status == KIN_SUCCESS, ExcMessage("Error initializing KINSOL. KINInit failed."));

  // set the level of verbosity
  status = KINSetPrintLevel(kin_mem, verbosity);
  Assert(status == KIN_SUCCESS, ExcMessage("Error intializing KINSOL. KINSetPrintLevel failed."));

  status = KINSetMaxSetupCalls(kin_mem, mbset);
  Assert(status == KIN_SUCCESS, ExcMessage("Error intializing KINSOL. KINSetMaxSetupCalls failed."));


  if (use_internal_solver)
    status = KINDense(kin_mem, system_size );
  else
    {
      KINMem KIN_mem;
      KIN_mem = (KINMem) kin_mem;
      KIN_mem->kin_lsetup = kinsol_setup_jacobian<VEC>;
      KIN_mem->kin_lsolve = kinsol_solve_linear_system<VEC>;
      KIN_mem->kin_setupNonNull = true;
    }
  is_initialized = true;

  (void)status;
}

// Run the solver:
template <typename VEC>
int KINSOLInterface<VEC>::solve( VEC &sol )
{
  int status=0;
  // Check initialization:
  if (!is_initialized)
    {
      initialize_solver( sol );
    }
  // Check initialization of the scaling. If scaling is not initialized we modify the
  // call to spare allocating a vector.
  if ( !scaling_is_set )
    {
#ifdef DEAL_II_WITH_MPI
      u_scale    = N_VNew_Parallel(communicator, local_system_size, system_size);
      N_VConst_Parallel( 1.e0, u_scale );
      f_scale    = N_VNew_Parallel(communicator, local_system_size, system_size);
      N_VConst_Parallel( 1.e0, f_scale );

#else
      u_scale    = N_VNew_Serial(system_size);
      N_VConst_Serial( 1.e0, u_scale );
      f_scale = N_VNew_Serial(system_size);
      N_VConst_Serial( 1.e0, f_scale );
#endif
    }


  // call to KINSol:
  if (strategy == "newton")
    status = KINSol(kin_mem, this->solution, KIN_NONE, u_scale, f_scale);
  if (strategy == "global_newton")
    status = KINSol(kin_mem, this->solution, KIN_LINESEARCH, u_scale, f_scale);
  if (strategy == "fixed_point")
    status = KINSol(kin_mem, this->solution, KIN_FP, u_scale, f_scale);
  if (strategy == "picard")
    status = KINSol(kin_mem, this->solution, KIN_PICARD, u_scale, f_scale);

  if (status == KIN_MXNEWT_5X_EXCEEDED)
    {
      // we get here when the inequality
      // norm_L2(u_scale*newton_update) > 0.99*mxnewtstep
      // is satisfied for 5 consecutive steps.
      // Such a failure may mean that norm_L2(f_scale*F(u))
      // asymptotes from above to a positive value,
      // or the real scalar mxnewtstep is too small.
      // So we rescale the mxnewtstep to the residual norm
      // and give to kinsol another try
      auto res = create_new_vector();
      KINMem KIN_mem;
      KIN_mem = (KINMem) kin_mem;
      copy(*res, KIN_mem->kin_fval);
      KINSetMaxNewtonStep(kin_mem, res->l2_norm());

      pcout << "Don't worry. This might be only a problem of scaling... Let's try."
            << std::endl;

      if (strategy == "newton")
        status = KINSol(kin_mem, this->solution, KIN_NONE, u_scale, f_scale);
      if (strategy == "global_newton")
        status = KINSol(kin_mem, this->solution, KIN_LINESEARCH, u_scale, f_scale);

    }
  AssertThrow(status >= 0 , ExcMessage("KINSOL did not converge. You might try with a different strategy."));

  copy( sol, this->solution );
  return status;

}

template <typename VEC>
void KINSOLInterface<VEC>::set_scaling_vectors( const VEC &uscale, const VEC &fscale )
{

  // Check initialization:
  AssertThrow(is_initialized,
              ExcMessage("KINSOLInterface is not initialized when "
                         "calling set_scaling_vectors!\n"
                         "you need to call initialize_solver() first "
                         "in order to make kinsol aware of the dimension of the system."));
  // initialize the scaling vectors:
#ifdef DEAL_II_WITH_MPI
  u_scale    = N_VNew_Parallel(communicator, local_system_size, system_size);
  f_scale    = N_VNew_Parallel(communicator, local_system_size, system_size);
#else
  u_scale    = N_VNew_Serial(system_size);
  f_scale    = N_VNew_Serial(system_size);
#endif
  // copy:
  copy(u_scale, uscale );
  copy(f_scale, fscale );
  // scaling is now initialized:
  scaling_is_set = true;

}

template <typename VEC>
void KINSOLInterface<VEC>::set_constraint_vector( const VEC &constraint )
{

  int status;
  // Check initialization:
  AssertThrow(is_initialized,
              ExcMessage("KINSOLInterface is not initialized "
                         "when calling initialize_constraint_vector!"));
  // copy the vector to an N_Vector:
#ifdef DEAL_II_WITH_MPI
  N_Vector constraint_v    = N_VNew_Parallel(communicator, local_system_size, system_size);
#else
  N_Vector constraint_v    = N_VNew_Serial(system_size);
#endif
  copy(constraint_v, constraint );
  // call:
  status = KINSetConstraints( kin_mem, constraint_v);
  // check status:
  Assert(status == KIN_SUCCESS, ExcMessage("Error initializing KINSOL. KINSetConstraints in initialize_constraint_vector failed."));
  (void)status;

}

D2K_NAMESPACE_CLOSE

template class deal2lkit::KINSOLInterface<BlockVector<double> >;

#ifdef DEAL_II_WITH_MPI

#ifdef DEAL_II_WITH_TRILINOS
template class deal2lkit::KINSOLInterface<TrilinosWrappers::MPI::Vector>;
template class deal2lkit::KINSOLInterface<TrilinosWrappers::MPI::BlockVector>;
#endif

#ifdef DEAL_II_WITH_PETSC
template class deal2lkit::KINSOLInterface<PETScWrappers::MPI::Vector>;
template class deal2lkit::KINSOLInterface<PETScWrappers::MPI::BlockVector>;
#endif

#endif

#endif
