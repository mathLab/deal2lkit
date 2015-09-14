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
#include <deal2lkit/fe_values_cache.h>
#include <deal2lkit/parsed_grid_generator.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/vector_tools.h>

int main ()
{
  initlog();

  ParsedGridGenerator<2> gg;

  ParameterAcceptor::initialize();

  auto tria = gg.serial();

  // Exact for quadratic functions
  FESystem<2> fe(FE_Q<2>(2), 2);

  DoFHandler<2> dh(*tria);

  tria->refine_global(1);
  dh.distribute_dofs(fe);

  Vector<double> sol(dh.n_dofs());
  Vector<double> sol_dot(dh.n_dofs());

  // x*y
  Tensor<1,2> e;
  e[0] = 1;
  e[1] = 1;

  Functions::Monomial<2> f(e, 2);

  VectorTools::interpolate(dh, f, sol);

  FEValuesCache<2,2> cache(StaticMappingQ1<2>::mapping,
                           fe, QGauss<2>(3),
                           update_values|
                           update_gradients|
                           update_quadrature_points|
                           update_JxW_values,
                           QGauss<1>(1),
                           update_default);
  double dummy = 1.0;


  std::vector<double> error(5, 0.0);

  FEValuesExtractors::Scalar p(0);
  FEValuesExtractors::Vector u(0);

  for (auto cell : dh.active_cell_iterators())
    {
      cache.reinit(cell);
      cache.cache_local_solution_vector("solution", sol, dummy);

      auto &us = cache.get_values("solution", "u", u, dummy);
      auto &ps = cache.get_values("solution", "p", p, dummy);

      auto &div_us = cache.get_divergences("solution", "u", u, dummy);
      auto &grad_us = cache.get_gradients("solution", "u", u, dummy);
      auto &grad_ps = cache.get_gradients("solution", "p", p, dummy);

      auto &Fs = cache.get_deformation_gradients("solution", "u", u, dummy);
      auto &sym_grad_us = cache.get_symmetric_gradients("solution", "u", u, dummy);

      auto &p = cache.get_quadrature_points();
      auto &JxW = cache.get_JxW_values();

      for (unsigned int q=0; q<us.size(); ++q)
        {
          for (unsigned int d=0; d<2; ++d)
            error[0] += (us[q][d]-f.value(p[q], d))*JxW[q];
          error[1] += (ps[q] - f.value(p[q], 0))*JxW[q];
          // div(u) = y+x
          error[2] += (div_us[q] - p[q][1] - p[q][0])*JxW[q];
          error[3] += (grad_ps[q][0]-p[q][1])*JxW[q];
          error[4] += (grad_ps[q][1]-p[q][0])*JxW[q];
        }
    }

  for (unsigned int i=0; i<error.size(); ++i)
    deallog << "Error[" << i << "]: "<<  error[i] << std::endl;

}
