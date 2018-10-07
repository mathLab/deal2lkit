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

// Test error handler with two components, one name.
// AddUp implicitly set.

#include <deal.II/base/function_lib.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal2lkit/error_handler.h>
#include <deal2lkit/parsed_grid_generator.h>
#include <deal2lkit/utilities.h>

#include "../tests.h"


using namespace deal2lkit;

int main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();
#else
  initlog();
#endif
  deallog.depth_console(0);

  ParsedGridGenerator<2, 3> gg;
  ErrorHandler<>            eh("", "u,u", "L2;L2"); // Only one table

  ParameterAcceptor::initialize();

#ifdef DEAL_II_WITH_MPI
#  include "mpi.h"
  auto tria = SP(gg.distributed(MPI_COMM_WORLD));
#else
  auto tria = gg.serial();
#endif

  FESystem<2, 3>   fe(FE_Q<2, 3>(1), 2);
  DoFHandler<2, 3> dh(*tria);

  for (unsigned int i = 0; i < 5; ++i)
    {
      tria->refine_global(1);
      dh.distribute_dofs(fe);
      Vector<double> sol(dh.n_dofs());
      VectorTools::interpolate(dh, Functions::CosineFunction<3>(2), sol);
      eh.error_from_exact(dh, sol, Functions::CosineFunction<3>(2));
    }
#ifdef DEAL_II_WITH_MPI
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
#endif
      eh.output_table(deallog.get_file_stream());
#ifdef DEAL_II_WITH_MPI
    }
#endif
}
