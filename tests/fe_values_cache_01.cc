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
#include "fe_values_cache.h"
#include "parsed_grid_generator.h"
#include <deal.II/base/function_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/vector_tools.h>

int main ()
{
  initlog();

  ParsedGridGenerator<2> gg;

  ParameterAcceptor::initialize();

  auto tria = gg.serial();

  FE_Q<2> fe(1);
  DoFHandler<2> dh(*tria);

  tria->refine_global(1);
  dh.distribute_dofs(fe);

  Vector<double> sol(dh.n_dofs());
  Vector<double> sol_dot(dh.n_dofs());

  Functions::CosineFunction<2> f(1);
  Functions::ExpFunction<2> g;

  VectorTools::interpolate(dh, f, sol);
  VectorTools::interpolate(dh, g, sol_dot);

  FEValuesCache<2,2> cache(StaticMappingQ1<2>::mapping,
                           fe, QGauss<2>(3),
                           update_values|
                           update_quadrature_points|
                           update_JxW_values,
                           QGauss<1>(1),
                           update_default);
  double dummy = 1.0;

  double error0 = 0;
  double error1 = 0;

  for (auto cell : dh.active_cell_iterators())
    {
      cache.reinit(cell);
      cache.cache_local_solution_vector("solution", sol, dummy);
      cache.cache_local_solution_vector("solution_dot", sol_dot, dummy);

      auto &us = cache.get_values("solution", dummy);
      auto &us_dot = cache.get_values("solution_dot", dummy);

      auto &p = cache.get_quadrature_points();
      auto &JxW = cache.get_JxW_values();

      for (unsigned int q=0; q<us.size(); ++q)
        {
          double diff0 = f.value(p[q]) - us[q][0];
          error0 += diff0*diff0*JxW[q];

          double diff1 = g.value(p[q]) - us_dot[q][0];
          error1 += diff1*diff1*JxW[q];

        }
    }

  deallog << "Error cosine: " << error0 << std::endl;
  deallog << "Error exp   : " << error1 << std::endl;

}
