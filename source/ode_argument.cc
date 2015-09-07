#ifdef DEAL_II_SAK_WITH_SUNDIALS

#include "ode_argument.h"
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#ifdef DEAL_II_WITH_TRILINOS
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#endif


#include "utilities.h"

template<typename VEC>
bool OdeArgument<VEC>::solver_should_restart(const double,
                                             const VEC &,
                                             const VEC &,
                                             const unsigned int,
                                             const double)
{
  return false;
}


template<typename VEC>
VEC &OdeArgument<VEC>::differential_components() const
{
  static shared_ptr<VEC> tmp = create_new_vector();
  static bool init = true;
  if (init == true)
    {
      tmp->add(1.0);
      init = false;
    }
  return (*tmp);
}

template<typename VEC>
VEC &OdeArgument<VEC>::get_local_tolerances() const
{
  static shared_ptr<VEC> tmp = create_new_vector();
  static bool init = true;
  if (init == true)
    {
      tmp->add(1.0);
      init = false;
    }
  return *tmp;
}

template class OdeArgument<Vector<double> >;
template class OdeArgument<BlockVector<double> >;


#ifdef DEAL_II_WITH_TRILINOS
template class OdeArgument<TrilinosWrappers::MPI::BlockVector>;
#endif

#endif
