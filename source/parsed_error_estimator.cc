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
                      const std::string &default_name,
                      const std::string &_estimator_strategy)
  :
  ParameterAcceptor(name),
  _n_components(n_components),
  mask(true),
  estimator_strategy(_estimator_strategy)
{
  parse_parameters_call_back();
}

template <int dim, int spacedim, int n_components>
void
ParsedErrorEstimator<dim, spacedim, n_components>::
declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &estimator_strategy,
                "Estimator strategy", estimator_strategy,
                Patterns::Selection("cell_diameter_over_24|face_diameter_over_twice_max_degree"),
                "- cell_diameter_over_24 : Kelly error estimator with the factor h/24.\n"
                "- face_diameter_over_twice_max_degree : the boundary residual estimator with the factor h_F/(2max(p+,p−)).");
}

template <int dim, int spacedim, int n_components>
template<typename InputVector>
void
ParsedErrorEstimator<dim, spacedim, n_components>::
compute_estimator(
  const DoFHandler<dim,spacedim> &dof,
  const Quadrature< dim-1 > &   quadrature,
  const InputVector    &solution,
  const typename FunctionMap< spacedim >::type   &neumann_bc,
  const Function< spacedim >  *coefficients,
  const unsigned int  n_threads,
  const types::subdomain_id subdomain_id,
  const types::material_id material_id
)
{
  KellyErrorEstimator< dim, spacedim >::
  estimate  (
    dof,
    quadrature,
    neumann_bc,
    solution,
    estimated_error_per_cell,
    ComponentMask(mask),
    coefficients,
    n_threads,
    subdomain_id,
    material_id,
    estimator_strategy
  );
};

template <int dim, int spacedim, int n_components>
double
ParsedErrorEstimator<dim, spacedim, n_components>::
linfty_norm()
{
  return estimated_error_per_cell.linfty_norm();
};

template <int dim, int spacedim, int n_components>
void
ParsedErrorEstimator<dim, spacedim, n_components>::
set_mask(const std::vector<bool> new_mask)
{
  mask = new_mask;
};


D2K_NAMESPACE_CLOSE
