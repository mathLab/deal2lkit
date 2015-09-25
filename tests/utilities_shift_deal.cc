#include "tests.h"
#include <deal2lkit/utilities.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>
#ifdef DEAL_II_WITH_MPI
#include <deal.II/lac/parallel_vector.h>
#endif


int main (int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);
#endif

#ifdef DEAL_II_WITH_MPI
  Assert(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 2, ExcNotImplemented());
  IndexSet index1;
  index1.set_size(10);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    for (unsigned int i = 0; i<5; ++i)
      index1.add_index(i);
  else
    for (unsigned int i = 5; i<10; ++i)
      index1.add_index(i);
  index1.compress();
  parallel::distributed::Vector<double> v1(index1, MPI_COMM_WORLD);
  for (IndexSet::size_type i=0; i < index1.n_elements(); ++i)
    v1[index1.nth_index_in_set(i)] = 1.;

  vector_shift(v1,1.);
  v1.print(std::cout);
#endif
}
