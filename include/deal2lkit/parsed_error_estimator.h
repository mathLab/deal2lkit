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
template <int dim, int spacedim=dim, int n_components=1>
class ParsedErrorEstimator : public ParameterAcceptor
{
public:
  /**
   * TODO:
   */
  ParsedErrorEstimator( const std::string &name =
                              "",
                        const std::string &estimator =
                              "kelly",
                        const std::string &string_component_mask =
                              "1,1,0",
                        const std::string &string_subdomain_id =
                              "numbers::invalid_subdomain_id",
                        const std::string &string_material_id =
                              "numbers::invalid_material_id",
                        const std::string &string_estimator_strategy =
                              "cell_diameter_over_24",
                        const string_n_threads =
                              "0");

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
    const DoFHandler<dim,spacedim>  &dof,
    const InputVector               &solution,
    const Vector<float>             &estimated_error_per_cell;
  );

private:

  std::string string_estimator_strategy;
  std::string string_subdomain_id;
  std::string string_material_id;
  std::string string_component_mask;
  std::string string_n_threads;

  KellyErrorEstimator::Strategy   estimator_strategy;
  unsigned int                    subdomain_id;
  unsigned int                    material_id;
  unsigned int                    n_threads;
  std::vector<bool>               mask;

};

D2K_NAMESPACE_CLOSE

#endif
