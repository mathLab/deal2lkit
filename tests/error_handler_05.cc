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

#include <deal2lkit/error_handler.h>
#include <deal2lkit/parsed_grid_generator.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>

int main ()
{
  initlog();

  ParsedGridGenerator<2,3> gg;
  ErrorHandler<> eh("Error", "u","L2, Linfty, H1, Custom"); // Only one table

  ParameterAcceptor::initialize();

  auto tria = gg.serial();

  FE_Q<2,3> fe(1);
  DoFHandler<2,3> dh(*tria);

  for (unsigned int i=0; i<5; ++i)
    {
      tria->refine_global(1);
      dh.distribute_dofs(fe);
      Vector<double> sol(dh.n_dofs());
      VectorTools::interpolate(dh, Functions::CosineFunction<3>(1), sol);
      eh.error_from_exact(dh, sol, Functions::CosineFunction<3>(1));
      eh.custom_error(  [&dh](unsigned int)
      {
        return 1.0/dh.get_tria().n_active_cells();
      }, dh);
    }
  eh.output_table(deallog.get_file_stream());
}
