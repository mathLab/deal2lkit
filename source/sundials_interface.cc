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


#include <deal2lkit/sundials_interface.h>
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


D2K_NAMESPACE_OPEN

template<typename VEC>
double SundialsInterface<VEC>::vector_norm(const VEC &vector) const
{
  return vector.l2_norm();
}

template<typename VEC>
bool SundialsInterface<VEC>::solver_should_restart(const double,
                                                   const unsigned int,
                                                   const double,
                                                   VEC &,
                                                   VEC &)
{
  return false;
}


template<typename VEC>
VEC &SundialsInterface<VEC>::differential_components() const
{
  static shared_ptr<VEC> tmp = create_new_vector();
  static bool init = true;
  if (init == true)
    {
      vector_shift(*tmp, 1.0);
      init = false;
    }
  return (*tmp);
}

template<typename VEC>
VEC &SundialsInterface<VEC>::get_local_tolerances() const
{
  static shared_ptr<VEC> tmp = create_new_vector();
  static bool init = true;
  if (init == true)
    {
      vector_shift(*tmp, 1.0);
      init = false;
    }
  return *tmp;
}

template<typename VEC>
int SundialsInterface<VEC>::jacobian_vmult(const VEC &, VEC &) const
{
  Assert(false, ExcPureFunctionCalled());
  int i=0;
  return i;
}


template<typename VEC>
void SundialsInterface<VEC>::get_lumped_mass_matrix(VEC &diag) const
{
  diag=1;
}

D2K_NAMESPACE_CLOSE

template class deal2lkit::SundialsInterface<Vector<double> >;
template class deal2lkit::SundialsInterface<BlockVector<double> >;

#ifdef DEAL_II_WITH_MPI
#ifdef DEAL_II_WITH_PETSC
template class deal2lkit::SundialsInterface<PETScWrappers::MPI::Vector>;
template class deal2lkit::SundialsInterface<PETScWrappers::MPI::BlockVector>;
#endif

#ifdef DEAL_II_WITH_TRILINOS
template class deal2lkit::SundialsInterface<TrilinosWrappers::MPI::Vector>;
template class deal2lkit::SundialsInterface<TrilinosWrappers::MPI::BlockVector>;
#endif
#endif

#endif
