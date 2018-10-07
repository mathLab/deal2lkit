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

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal2lkit/error_handler.h>
#include <deal2lkit/parsed_grid_generator.h>

#include "../tests.h"


using namespace deal2lkit;

int main()
{
  initlog();

  ParsedGridGenerator<2, 3> gg;
  ErrorHandler<>            eh; // Only one table

  ParameterAcceptor::initialize();

  auto tria = gg.serial();

  FE_Q<2, 3>       fe(1);
  DoFHandler<2, 3> dh(*tria);

  for (unsigned int i = 0; i < 5; ++i)
    {
      tria->refine_global(1);
      dh.distribute_dofs(fe);
      Vector<double> sol(dh.n_dofs());
      VectorTools::interpolate(dh, Functions::CosineFunction<3>(1), sol);
      eh.error_from_exact(dh, sol, Functions::CosineFunction<3>(1));
    }
  eh.output_table(deallog.get_file_stream());
}
