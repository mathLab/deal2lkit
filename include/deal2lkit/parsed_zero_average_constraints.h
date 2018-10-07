//-----------------------------------------------------------
//
//    Copyright (C) 2016 by the deal2lkit authors
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

#ifndef _d2k_parsed_zero_average_constraints_h
#define _d2k_parsed_zero_average_constraints_h

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parsed_function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>

#include <deal.II/lac/constraint_matrix.h>

#include <deal2lkit/config.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/utilities.h>

#include <algorithm>
#include <map>


using namespace dealii;


D2K_NAMESPACE_OPEN

/**
 * ParsedZeroAverageConstraints class. It allows you to set a zero
 * average constraint to a ConstraintMatrix. The variable to be
 * constrained (e.g. the pressure for Stokes problem) can be specified
 * in the parameter file either by its component name (e.g. "p") or by
 * the component number.
 *
 * The zero mean value can be set on the whole domain, or on the
 * boundary.
 *
 * A typical usage of this class is as follows:
 *
 * @code
 * ConstraintMatrix cm;
 * ParsedZeroAverageConstraints<dim> pnac("Parsed Zero Average Constraints",
 *              dim+1,
 *            (dim==2?"u,u,p":"u,u,u,p"),
 *            "p");
 *
 *
 *
 * ParameterAcceptor::initialize();
 * pnac.apply_zero_average_constraints(dof_handler,cm);
 * @endcode
 */



template <int dim, int spacedim = dim>
class ParsedZeroAverageConstraints : public ParameterAcceptor
{
public:
  /**
   * Constructor.
   *
   * It takes:
   * - the name for the section of the Parameter Handler to use
   *
   * - the number of components that this function has
   *
   * - the default names of known components that can be used instead
   *   of component numbers. If empty, in the parameter file will be
   *   filled with n_components 'u'
   *
   * - a list of components, separated by a "," (i.e. 0,p,5), that are
   *   to be constrained on the whole domain
   *
   * - a list of components, separated by a "," (i.e. 0,p,5), that are
   *   to be constrained on the boundary
   */
  ParsedZeroAverageConstraints(
    const std::string & name,
    const unsigned int &n_components                = 1,
    const std::string & component_names             = "",
    const std::string & default_components          = "",
    const std::string & default_boundary_components = "");


  /**
   * Compute the zero average constraints and apply them on the given
   * constraint matrix
   */
  void
  apply_zero_average_constraints(const DoFHandler<dim, spacedim> &dof_handler,
                                 ConstraintMatrix &constraints) const;


  /**
   * return the ComponentMask at boundary
   */
  ComponentMask get_boundary_mask() const;


  /**
   * return the ComponentMask
   */
  ComponentMask get_mask() const;


  /**
   * declare_parameters is inherithed by ParameterAcceptor
   */
  virtual void declare_parameters(ParameterHandler &prm);

  /**
   * parse_parameters_call_back is inherithed by ParameterAcceptor
   */
  virtual void parse_parameters_call_back();



  /// Wrong number of component mask
  DeclException2(ExcWrongComponent,
                 unsigned int,
                 unsigned int,
                 << "Wrong component number has been used: " << arg1
                 << " is not in the range [0, " << arg2 << ").");

  /// Wrong variable name
  DeclException2(ExcWrongVariable,
                 std::string,
                 std::vector<std::string>,
                 << "Wrong variabile name has been used: " << arg1
                 << " does not belong to the knwon variables: "
                 << print(unique(arg2)) << ".");

protected:
  void internal_zero_average_constraints(
    const DoFHandler<dim, spacedim> &dof_handler,
    const ComponentMask              mask,
    const bool                       at_boundary,
    ConstraintMatrix &               constraints) const;

  std::string              name;
  std::string              str_components;
  std::string              str_boundary_components;
  std::string              str_component_names;
  std::vector<std::string> _component_names;
  std::vector<std::string> boundary_components;
  std::vector<std::string> components;

  std::vector<bool> mask;
  std::vector<bool> boundary_mask;

  const unsigned int n_components;
};



D2K_NAMESPACE_CLOSE

#endif
