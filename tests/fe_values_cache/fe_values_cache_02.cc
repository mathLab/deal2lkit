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

#include <deal.II/numerics/vector_tools.h>

#include <deal2lkit/fe_values_cache.h>
#include <deal2lkit/parsed_grid_generator.h>

#include "../tests.h"

using namespace deal2lkit;


int
main()
{
  initlog();

  ParsedGridGenerator<2> gg;

  ParameterAcceptor::initialize();

  auto tria = gg.serial();

  FE_Q<2>       fe(1);
  DoFHandler<2> dh(*tria);

  tria->refine_global(1);
  dh.distribute_dofs(fe);

  Vector<double> sol(dh.n_dofs());
  Vector<double> sol_dot(dh.n_dofs());

  Functions::CosineFunction<2> f(1);
  Functions::ExpFunction<2>    g;

  VectorTools::interpolate(dh, f, sol);
  VectorTools::interpolate(dh, g, sol_dot);

  FEValuesCache<2, 2> cache(StaticMappingQ1<2>::mapping,
                            fe,
                            QGauss<2>(3),
                            update_values | update_quadrature_points |
                              update_JxW_values,
                            QGauss<1>(1),
                            update_default);
  double              dummy = 1.0;

  double error0 = 0;
  double error1 = 0;

  FEValuesExtractors::Scalar u(0);

  for (auto cell : dh.active_cell_iterators())
    {
      cache.reinit(cell);
      cache.cache_local_solution_vector("solution", sol, dummy);
      cache.cache_local_solution_vector("solution_dot", sol_dot, dummy);

      auto &us     = cache.get_values("solution", "u", u, dummy);
      auto &us_dot = cache.get_values("solution_dot", "p", u, dummy);

      auto &p   = cache.get_quadrature_points();
      auto &JxW = cache.get_JxW_values();

      for (unsigned int q = 0; q < us.size(); ++q)
        {
          double diff0 = f.value(p[q]) - us[q];
          error0 += diff0 * diff0 * JxW[q];

          double diff1 = g.value(p[q]) - us_dot[q];
          error1 += diff1 * diff1 * JxW[q];
        }
    }

  deallog << "Error cosine: " << error0 << std::endl;
  deallog << "Error exp   : " << error1 << std::endl;
}
