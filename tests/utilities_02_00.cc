//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

// Test if the "create directory" utilities function works as expected.

#include "tests.h"
#include <deal2lkit/utilities.h>
#include <deal.II/base/utilities.h>

#ifdef DEAL_II_WITH_MPI
#include <deal.II/base/mpi.h>
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
  if ( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      // std::string = dealii::Utilities::int_to_string(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
#endif
      std::system("rm -rf dsfas*");
      deallog << "--> " << create_directory("dsfas000") << std::endl;
      deallog << "--> " << create_directory("dsfas001") << std::endl;
      deallog << "--> " << create_directory("dsfas002") << std::endl;
      deallog << "--> " << dir_exists("dsfas002") << std::endl;
      // deallog << "--> " << get_next_available_directory_name("dsfas",3) << std::endl;
      // deallog << "--> " << create_directory(get_next_available_directory_name("dsfas",3)) << std::endl;
      // deallog << "--> " << create_directory("dsfas003") << std::endl;
      // deallog << "--> " << get_next_available_directory_name("dsfas",3) << std::endl;
      return 0;
#ifdef DEAL_II_WITH_MPI
    }
  else
    {
      return 0;
    }
#endif
}
