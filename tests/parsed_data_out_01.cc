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

#include "tests.h"
#include "parsed_data_out.h"

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
