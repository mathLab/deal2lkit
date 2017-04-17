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

#include <deal2lkit/parsed_error_estimator.h>

// #include <deal.II/base/logstream.h>
// #include <deal.II/base/quadrature_lib.h>
// #include <deal.II/base/utilities.h>
// #include <deal.II/base/conditional_ostream.h>
//
// #include <deal.II/grid/grid_tools.h>
//
// #include <deal.II/numerics/vector_tools.h>
// #include <deal.II/numerics/matrix_tools.h>
// #include <deal.II/numerics/data_out.h>
// #include <deal.II/fe/mapping_q.h>
// #include <deal.II/fe/fe.h>
// #include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>


D2K_NAMESPACE_OPEN

template <int dim, int spacedim, int n_components>
ParsedErrorEstimator<dim, spacedim, n_components>::
ParsedErrorEstimator( const std::string &name,
                      const std::string &estimator,
                      const std::string &string_subdomain_id,
                      const std::string &string_material_id,
                      const std::string &string_estimator_strategy,
                      const string_n_threads);
  :
  ParameterAcceptor(name),
  _n_components(n_components),
  mask(true),
  string_estimator_strategy(string_estimator_strategy)
{
  parse_parameters_call_back();
}

template <int dim, int spacedim, int n_components>
void
ParsedErrorEstimator<dim, spacedim, n_components>::
declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &string_estimator_strategy,
                "Estimator strategy", string_estimator_strategy,
                Patterns::Selection("cell_diameter_over_24|face_diameter_over_twice_max_degree"),
                "- cell_diameter_over_24 : Kelly error estimator with the factor h/24.\n"
                "- face_diameter_over_twice_max_degree : the boundary residual estimator with the factor h_F/(2max(p+,p−)).");
}

template <int dim, int spacedim, int n_components>
void
ParsedErrorEstimator<dim, spacedim, n_components>::
parse_parameters_call_back()
{
  if(string_estimator_strategy=="cell_diameter_over_24")
    estimator_strategy = cell_diameter_over_24;
  else
    estimator_strategy = face_diameter_over_twice_max_degree;
}

template <int dim, int spacedim, int n_components>
template<typename InputVector>
void
ParsedErrorEstimator<dim, spacedim, n_components>::
compute_estimator(
  const DoFHandler<dim,spacedim> &dof,
  const InputVector              &solution,
  const Vector<float>            &estimated_error_per_cell;
)
{
  KellyErrorEstimator<dim>::estimate (dh,
                                      QGauss<dim-1>(dh.get_fe().degree + 1),
                                      typename FunctionMap<dim>::type(),
                                      solution,
                                      estimated_error_per_cell,
                                      ComponentMask(mask),
                                      0,
                                      n_threads,
                                      subdomain_id,
                                      material_id,
                                      strategy);
};

D2K_NAMESPACE_CLOSE
