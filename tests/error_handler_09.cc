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

// Test error handler with two components, one name.
// AddUp implicitly set.

#include "tests.h"
#include "utilities.h"
#include "error_handler.h"
#include "parsed_grid_generator.h"
#include <deal.II/base/function_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>

int main (int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();
#else
  initlog();
#endif
  deallog.depth_console (0);

  ParsedGridGenerator<2,3> gg;
  ErrorHandler<> eh("", "u,u", "L2;L2"); // Only one table

  ParameterAcceptor::initialize();

#ifdef DEAL_II_WITH_MPI
#include "mpi.h"
  auto tria = SP(gg.distributed(MPI_COMM_WORLD));
#else
  auto tria = gg.serial();
#endif

  FESystem<2,3> fe(FE_Q<2,3>(1), 2);
  DoFHandler<2,3> dh(*tria);

  for (unsigned int i=0; i<5; ++i)
    {
      tria->refine_global(1);
      dh.distribute_dofs(fe);
      Vector<double> sol(dh.n_dofs());
      VectorTools::interpolate(dh, Functions::CosineFunction<3>(2), sol);
      eh.error_from_exact(dh, sol, Functions::CosineFunction<3>(2));
    }
#ifdef DEAL_II_WITH_MPI
  if ( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
#endif
      eh.output_table(deallog.get_file_stream());
#ifdef DEAL_II_WITH_MPI
    }
#endif
}
