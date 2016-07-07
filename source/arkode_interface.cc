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

#include <deal2lkit/arkode_interface.h>

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


D2K_NAMESPACE_OPEN


// anonymus namespace to hide helper functions required
// to interface ARKode of Sundials
namespace
{
  template<typename VEC>
  int arkode_explicit_rhs (double t,
                           N_Vector y,
                           N_Vector ydot,
                           void *user_data)
  {
    ARKodeInterface<VEC> &solver =
      *static_cast<ARKodeInterface<VEC> >(user_data);

    shared_ptr<VEC> loc_y = solver.create_new_vector();
    shared_ptr<VEC> loc_ydot = solver.create_new_vector();

    copy(*loc_y,y);
    copy(*loc_ydoy,ydot);

    int ret = solver.explicit_rhs(t,*loc_y,*loc_ydot);

    copy(ydot,loc_ydot);

    return ret;
  }

  template<typename VEC>
  int arkode_implicit_rhs (double t,
                           N_Vector y,
                           N_Vector ydot,
                           void *user_data)
  {
    ARKodeInterface<VEC> &solver =
      *static_cast<ARKodeInterface<VEC> >(user_data);

    shared_ptr<VEC> loc_y = solver.create_new_vector();
    shared_ptr<VEC> loc_ydot = solver.create_new_vector();

    copy(*loc_y,y);
    copy(*loc_ydoy,ydot);

    int ret = solver.implicit_rhs(t,*loc_y,*loc_ydot);

    copy(ydot,*loc_ydot);

    return ret;
  }

  template <typename VEC>
  int arkode_lsetup(ARKodeMem arkode_mem,
                    int /*convfail*/,
                    N_Vector ypred,
                    N_Vector /*fpred*/,
                    bool *jcurPtr,
                    N_Vector /*vtemp1*/,
                    N_Vector /*vtemp2*/,
                    N_Vector /*vtemp3*/)
  {
    ARKodeInterface<VEC> &solver =
      *static_cast<ARKodeInterface<VEC> >(arkode_mem->ark_user_data);

    shared_ptr<VEC> loc_y = solver.create_new_vector();

    copy(*loc_y,ypred);

    const double time = arkode_mem->ark_tn;

    int ret = solver.setup_jacobian(time,*loc_y);

    *jcurPtr = true;

    return ret;
  }

  template <typename VEC>
  int arkode_msetup(ARKodeMem arkode_mem,
                    N_Vector /*vtemp1*/,
                    N_Vector /*vtemp2*/,
                    N_Vector /*vtemp3*/)
  {
    ARKodeInterface<VEC> &solver =
      *static_cast<ARKodeInterface<VEC> >(arkode_mem->ark_user_data);

    const double time = arkode_mem->ark_tn;

    int ret = solver.mass_matrix(time);

    return ret;
  }

  // solve (M-gamma*J)=b
  template <typename VEC>
  int arkode_lsolve (ARKodeMem arkode_mem,
                     N_Vector b,
                     N_Vector /*weight*/,
                     N_Vector /*ycur*/,
                     N_Vector fcur)
  {
    ARKodeInterface<VEC> &solver =
      *static_cast<ARKodeInterface<VEC> >(arkode_mem->ark_user_data);

    shared_ptr<VEC> src = solver.create_new_vector();
    shared_ptr<VEC> dst = solver.create_new_vector();

    copy(*src,fcur);

    const double gamma = arkode_mem->ark_gamma;

    int ret = solver.solve_linear_system(gamma,*src,dst);

    copy(b,*dst);

    return ret;
  }

  // solve Mx=b, where M is the mass matrix
  template <typename VEC>
  int arkode_msolve (ARKodeMem arkode_mem,
                     N_Vector b,
                     N_Vector /*weight*/)
  {
    ARKodeInterface<VEC> &solver =
      *static_cast<ARKodeInterface<VEC> >(arkode_mem->ark_user_data);

    shared_ptr<VEC> src = solver.create_new_vector();
    shared_ptr<VEC> dst = solver.create_new_vector();

    copy(*src,b);

    int ret = solver.solve_mass_system(*src,dst);

    copy(b,*dst);

    return ret;
  }

} // close anonymous namespace



template <typename VEC>
ARKodeInterface<VEC>::ARKodeInterface(const std::string name,
                                      const MPI_Comm mpi_comm):
  ParameterAcceptor(name),
  communicator(Utilities::MPI::duplicate_communicator(mpi_comm)),
  pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_comm)==0),
  ark_mem(nullptr)
{
  set_functions_to_trigger_an_assert();
}

#else

template <typename VEC>
ARKodeInterface<VEC>::ARKodeInterface(const std::string name):
  ParameterAcceptor(name),
  is_initialized(false),
  pcout(std::cout),
  arkode_mem(nullptr)
{
  set_functions_to_trigger_an_assert();
}

#endif

template <typename VEC>
void
ARKodeInterface<VEC>::initialize(const double &t0,
                                 const VEC &y0)
{
  int status;
  if (ark_mem)
    {
      ARKodeFree(&ark_mem);
    }
  ark_mem = ARKodeCreate();
  system_size = y0.size();

#ifdef DEAL_II_WITH_MPI
  IndexSet is = y0.locally_owned_elements();
  local_system_size = is.n_elements();
  internal_solution    = N_VNew_Parallel(communicator,
                                         local_system_size,
                                         system_size);
#else
  internal_solution    = N_VNew_Serial(system_size);
#endif

  copy(internal_solution,y0);

  status = ARKodeSetUserData(ark_mem, (void *) this);

  status = ARKodeInit(ark_mem,
                      arkode_explicit_rhs,
                      arkode_implicit_rhs,
                      t0,
                      y0);

  status = ARKodeSStolerances(ark_mem, rtol, atol);

  ARKodeMem ARK_mem;
  ARK_mem = (ARKodeMem)ark_mem;

  ARK_mem->ark_lsetup = arkode_lsetup<VEC>;
  ARK_mem->ark_msetup = arkode_msetup<VEC>;

  ARK_mem->ark_lsolve = arkode_lsolve<VEC>;
  ARK_mem->ark_msolve = arkode_msolve<VEC>;

  ARK_mem->ark_setupNonNull = true;
  ARK_mem->ark_MassSetupNonNull = true;

}

template <typename VEC>
void
ARKodeInterface<VEC>::solve(VEC &solution)
{

}

template <typename VEC>
ARKodeInterface<VEC>::~ARKodeInterface()
{
  if (ark_mem)
    ARKodeFree(&ark_mem);
}

template <typename VEC>
void
ARKodeInterface<VEC>::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm,&rtol,"Relative tolerance", "1e-3",Patterns::Double(0));

  add_parameter(prm,&atol,"Absolute tolerance", "1e-4",Patterns::Double(0));

}

D2K_NAMESPACE_CLOSE

template class deal2lkit::ARKodeInterface<BlockVector<double> >;

#ifdef DEAL_II_WITH_MPI

#ifdef DEAL_II_WITH_TRILINOS
template class deal2lkit::ARKodeInterface<TrilinosWrappers::MPI::Vector>;
template class deal2lkit::ARKodeInterface<TrilinosWrappers::MPI::BlockVector>;
#endif

#ifdef DEAL_II_WITH_PETSC
template class deal2lkit::ARKodeInterface<PETScWrappers::MPI::Vector>;
template class deal2lkit::ARKodeInterface<PETScWrappers::MPI::BlockVector>;
#endif

#endif

#endif
