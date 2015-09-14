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

// Test if the "copy files" utilities function works as expected

#include "tests.h"
#include <deal2lkit/utilities.h>
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
      std::system("touch test1.txt");
      std::system("touch test2.txt");
      std::system("rm -rf new_folder*");
      deallog << "copy   --> " << copy_file("test1.txt", "tmp") << std::endl;
      deallog << "exist  --> " << file_exists("test1.txt") << std::endl;
      deallog << "exist  --> " << file_exists("test2.txt") << std::endl;
      std::string folder = get_next_available_directory_name("new_folder",3);
      deallog << "create --> " << folder << " - "<< create_directory(folder) << std::endl;
      deallog << "copy   --> " << copy_files("test1.txt test2.txt", folder) << std::endl;
      deallog << "exist  --> " << file_exists(folder+"/"+"test1.txt") << std::endl;
      deallog << "exist  --> " << file_exists(folder+"/"+"test2.txt") << std::endl;
#ifdef DEAL_II_WITH_MPI
    }
#endif
}
