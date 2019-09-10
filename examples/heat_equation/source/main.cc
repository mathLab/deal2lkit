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

#include <deal.II/base/logstream.h>

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/utilities.h>

#include <fstream>

#include "heat_ida.h"

using namespace dealii;

int
main(int argc, char **argv)
{
  deallog.depth_console(0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);


  MPI_Comm comm = MPI_COMM_WORLD;


  Heat<2> solver(comm);

  dealii::ParameterAcceptor::initialize("../source/heat_ida.prm",
                                        "used_parameters.prm");

  solver.run();
}
