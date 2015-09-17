#include "tests.h"
#include <deal2lkit/utilities.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>
#ifdef DEAL_II_WITH_PETSC
#include <deal.II/lac/petsc_parallel_vector.h>
#endif


int main (int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();
#else
  initlog();
#endif

#ifdef DEAL_II_WITH_MPI
#ifdef DEAL_II_WITH_PETSC
  Assert(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 2, ExcNotImplemented());
  IndexSet index1;
  index1.set_size(10);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    for (unsigned int i = 0; i<5; ++i)
      index1.add_index(i);
  else
    for (unsigned int i = 5; i<10; ++i)
      index1.add_index(i);

  PETScWrappers::MPI::Vector v1(index1, MPI_COMM_WORLD);
  for (auto i : v1.locally_owned_elements())
    v1[i] = 1.;
  v1.compress(VectorOperation::insert);

  vector_shift(v1,1.);
  v1.compress(VectorOperation::add);
  Vector<double> foo(v1);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      for (auto i : foo.locally_owned_elements())
        deallog<<i<<" "<<foo[i]<<std::endl;
    }
#endif
#endif
}
