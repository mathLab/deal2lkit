#include "heat_ida.h"


int main(int argc, char **argv)
{
  //Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);


  MPI_Comm comm = MPI_COMM_WORLD;


  int numprocs = Utilities::MPI::n_mpi_processes(comm);
  int myid = Utilities::MPI::this_mpi_process(comm);
//
//  std::cout << "Process " << getpid() << " is " << myid
//            << " of " << numprocs << " processes" << std::endl;
//  if (myid == 0) system("read -p \"Press [Enter] key to start debug...\"");
//
//

  typedef TrilinosWrappers::MPI::Vector VEC;

  Heat<2> solver(comm);

  ParameterAcceptor::initialize("heat_ida.prm", "used_parameters.prm");

  solver.run();

  return (0);
}
