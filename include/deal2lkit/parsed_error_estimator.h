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

#ifndef _d2k_parsed_solver_h
#define _d2k_parsed_solver_h

#include <deal.II/base/numbers.h>

#include <deal2lkit/parameter_acceptor.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/function.h>
#include <deal.II/dofs/function_map.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/lac/vector.h>

using namespace dealii;
using namespace deal2lkit;

D2K_NAMESPACE_OPEN

/**
 * TODO:
 */
template <int dim, int spacedim=dim, int n_components=dim>
class ParsedErrorEstimator : public ParameterAcceptor
{
public:
  /**
   * TODO:
   */
  ParsedErrorEstimator( const std::string &name,
                        const std::string &default_name,
                        const std::string &_estimator_strategy =
                          "cell_diameter_over_24");

  /**
   * TODO:
   */
  virtual void declare_parameters(ParameterHandler &prm);

  /**
   * TODO:
   */
  virtual void parse_parameters_call_back();

  /**
   * TODO:
   */
  template<typename InputVector>
  void compute_estimator(
    const DoFHandler<dim,spacedim> &dof,
    const Quadrature< dim-1 > &   quadrature,
    const InputVector    &solution,
    const typename FunctionMap< spacedim >::type   &neumann_bc =
      typename FunctionMap<dim>::type(),
    const Function< spacedim >  *coefficients =
      0,
    const unsigned int  n_threads =
      numbers::invalid_unsigned_int,
    const types::subdomain_id   subdomain_id =
      numbers::invalid_subdomain_id,
    const types::material_id  material_id =
      numbers::invalid_material_id
  );

  /**
   * TODO:
   */
  double linfty_norm();

  /**
   * TODO:
   */
  void set_mask(const std::vector<bool> new_mask);

private:

  const unsigned int _n_components;
  std::string estimator_strategy;

  double estimator_linfty_norm;
  Vector<float> estimated_error_per_cell;
  std::vector<bool> mask;

};

D2K_NAMESPACE_CLOSE

#endif
