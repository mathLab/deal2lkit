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

// Test error handler with two components, two names.
// No AddUp explicitly set.
#include "tests.h"

#include "error_handler.h"
#include "parsed_grid_generator.h"
#include <deal.II/base/function_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_handler.h>

int main ()
{
  initlog();

  ParsedGridGenerator<2,3> gg;
  ErrorHandler<> eh("", "u,p", "L2;L2"); // Only one table

  ParameterAcceptor::initialize();

  auto tria = gg.serial();

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
  eh.output_table(deallog.get_file_stream());
}
