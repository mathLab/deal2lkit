
#include <deal2lkit/sundials_interface.h>
#ifdef D2K_WITH_SUNDIALS
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/generic_linear_algebra.h>

#ifdef DEAL_II_WITH_PETSC
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#endif

#ifdef DEAL_II_WITH_TRILINOS
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#endif


#include <deal2lkit/utilities.h>

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

template class SundialsInterface<Vector<double> >;
template class SundialsInterface<BlockVector<double> >;

#ifdef DEAL_II_WITH_PETSC
template class SundialsInterface<PETScWrappers::MPI::Vector>;
template class SundialsInterface<PETScWrappers::MPI::BlockVector>;
#endif

#ifdef DEAL_II_WITH_TRILINOS
template class SundialsInterface<TrilinosWrappers::MPI::Vector>;
template class SundialsInterface<TrilinosWrappers::MPI::BlockVector>;
#endif

#endif
