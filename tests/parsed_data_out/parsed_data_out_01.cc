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

#include "../tests.h"
#include <deal2lkit/parsed_data_out.h>


using namespace deal2lkit;

int main (int argc, char **argv)
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,1);
#endif
  initlog();
  ParsedDataOut<2,2> pp("Test", "gnuplot");

  ParameterAcceptor::initialize();

  ParameterAcceptor::prm.log_parameters(deallog);
}
