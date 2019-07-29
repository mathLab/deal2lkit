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

// Test if the "test utilities and ParsedDataOut" utilities function works as
// expected

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal2lkit/parsed_data_out.h>
#include <deal2lkit/utilities.h>

#include "../tests.h"


using namespace deal2lkit;

int
main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);
  mpi_initlog();
#else
  initlog();
#endif

  std::system("rm -rf solution*");

  ParsedDataOut<2, 2> pp(
    "test", "vtu", 1, "solution/run", "solution", "", MPI_COMM_WORLD);

  dealii::ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  std::system("touch solution/run000/test.txt");
  deallog << "exists ? = " << file_exists("solution/run000/test.txt")
          << std::endl;
}
