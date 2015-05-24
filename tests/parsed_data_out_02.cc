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

// Test if the "test utilities and ParsedDataOut" utilities function works as expected

#include "tests.h"
#include "utilities.h"
#include "parsed_data_out.h"
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

  std::system("rm -rf solution*");

  ParsedDataOut<2,2> pp("test",
                        "vtu",
                        "solution/run",
                        "solution",
                        MPI_COMM_WORLD);

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  std::system("touch solution/run000/test.txt");
  deallog << "exists ? = " << file_exists("solution/run000/test.txt") << std::endl;
}
