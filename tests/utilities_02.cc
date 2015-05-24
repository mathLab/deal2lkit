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

// Test if the "print" utilities function works as expected. This is
// somewhat the inverse of the dealii::Utilities::split_string_list(),
// but it works with arbitrary objects that support operator<<();

#include "tests.h"
#include "utilities.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>

int main (int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);
  mpi_initlog();
#else
  initlog();
#endif

#ifdef DEAL_II_WITH_MPI
  if ( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
#endif
      std::system("rm -rf dsfas*");
      deallog << "--> " << create_directory("dsfas000") << std::endl;
      deallog << "--> " << create_directory("dsfas001") << std::endl;
      deallog << "--> " << create_directory("dsfas002") << std::endl;
      deallog << "--> " << get_next_available_directory_name("dsfas",3) << std::endl;
      deallog << "--> " << create_directory(get_next_available_directory_name("dsfas",3)) << std::endl;
      deallog << "--> " << create_directory("dsfas004") << std::endl;
      deallog << "--> " << get_next_available_directory_name("dsfas",3) << std::endl;
#ifdef DEAL_II_WITH_MPI
    }
#endif
}
