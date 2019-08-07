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

#ifndef d2k_parsed_dirichlet_bcs_h
#define d2k_parsed_dirichlet_bcs_h

#include <deal.II/base/exceptions.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal2lkit/config.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_mapped_functions.h>


D2K_NAMESPACE_OPEN

/**
 * Parsed DirichletBCs.
 * This class allows to set dirichlet boundary conditions,
 * in the strong formulation, from a parameter file.
 *
 * This class is derived from ParsedMappedFunctions.
 *
 * The VectorTools::interpolate_boundary_values and
 * VectorTools::project_boundary_values functions of
 * the deal.II library have been wrapped.
 *
 *
 * A typical usage of this class is the following
 *
 *
 * @code
 * // the following variables are set just for the sake of completeness
 * unsigned int dim = 2;
 * unsigned int spacedim = 2;
 * unsigned int n_components = 3; // number of components of the problem
 *
 * // create parsed_dirichlet object
 *
 * ParsedDirichletBCs<dim,spacedim>
 *    parsed_dirichlet("Dirichlet BCs", // name for the section of the Parameter
 * Handler to use n_components,    // number of components of the problem
 *                     "u,u,p",         // names of known components that can be
 * used instead of component numbers "0=u % 1=2 % 3=u.N % 6=ALL", // boundary_id
 * = component;other_component % other_id = comp; other_comp "0=x;y;0 % 1=0;0;0
 * % 3=x;2;0 % 6=y*k;0;k", // boundary_id = expression % other_id =
 * other_expression "k=1"); // list of constants that can be used in the above
 * epressions
 * ...
 * QGauss<dim-1> quadrature(2);
 * MappingQ1<dim> mapping;
 * std::map<types::global_dof_index,double> dirichlet_dofs;
 * AffineConstraints<double>  constraints;
 * ...
 * // the following functions apply the boundary conditions for the
 * // boundary ids 0,1,6
 * parsed_dirichlet.interpolate_boundary_values(dof_handler,dirichlet_dofs);
 * ...
 * parsed_dirichlet.interpolate_boundary_values(dof_handler,constraints);
 * ...
 * parsed_dirichlet.interpolate_boundary_values(mapping,dof_handler,constraints);
 * ...
 * parsed_dirichlet.interpolate_boundary_values(dof_handler,quadrature,dirichlet_dofs);
 * ...
 * parsed_dirichlet.interpolate_boundary_values(dof_handler,quadrature,constraints);
 * ...
 * parsed_dirichlet.interpolate_boundary_values(mapping,dof_handler,quadrature,constraints);
 * ...
 *
 * // in order to apply the boundary conditions also to the
 * // boundary id 3, where the normal component is set,
 * // the following function must be called
 *
 * parsed_dirichlet.compute_nonzero_normal_flux(dof_handler,constraints);
 *
 * @endcode
 */

template <int dim, int spacedim = dim>
class ParsedDirichletBCs : public ParsedMappedFunctions<spacedim>
{
public:
  /**
   * Constructor.
   *
   * It takes:
   * - the name for the section of the Parameter Handler to use
   *
   * - the names of known components that can be used instead of component
   * numbers
   *
   * - a list of ids and components where the boundary conditions must be
   * applied, which is a string with the following pattern boundary_id =
   * component;other_component % other_id = comp; other_comp
   *
   * - a list of ids and exprssions defined over the ids
   * (if this string is left empty, homogeneous boundary conditions are imposed
   * on the above specified ids and components)
   *
   * - list of constants that can be used in the above epressions
   *
   */
  ParsedDirichletBCs(const std::string & name            = "Dirichlet BCs",
                     const unsigned int &n_components    = 1,
                     const std::string & component_names = "",
                     const std::string & default_id_components = "0=ALL",
                     const std::string & default_id_functions  = "",
                     const std::string & default_constants     = "");

  /**
   * these method calls the method of the Parent class
   */
  virtual void
  declare_parameters(dealii::ParameterHandler &prm);

  /**
   * these method calls the method of the Parent class
   */
  virtual void
  parse_parameters_call_back();

  /**
   * This function must be called in order to apply the boundary conditions
   * to the AffineConstraints<double> .
   * It relies on the VectorTools::interpolate_boundary_values functions of the
   * deal.II library
   */
  void
  interpolate_boundary_values(
    const dealii::DoFHandler<dim, spacedim> &dof_handler,
    dealii::AffineConstraints<double> &      constraints) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the AffineConstraints<double> .
   * It relies on the VectorTools::interpolate_boundary_values functions of the
   * deal.II library
   */
  void
  interpolate_boundary_values(
    const dealii::Mapping<dim, spacedim> &   mapping,
    const dealii::DoFHandler<dim, spacedim> &dof_handler,
    dealii::AffineConstraints<double> &      constraints) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the AffineConstraints<double> .
   * It relies on the VectorTools::interpolate_boundary_values functions of the
   * deal.II library
   */
  void
  interpolate_boundary_values(
    const dealii::DoFHandler<dim, spacedim> &          dof_handler,
    std::map<dealii::types::global_dof_index, double> &d_dofs) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the AffineConstraints<double> .
   * It relies on the VectorTools::interpolate_boundary_values functions of the
   * deal.II library
   */
  void
  interpolate_boundary_values(
    const dealii::Mapping<dim, spacedim> &             mapping,
    const dealii::DoFHandler<dim, spacedim> &          dof_handler,
    std::map<dealii::types::global_dof_index, double> &d_dofs) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the AffineConstraints<double> .
   * It relies on the VectorTools::project_boundary_values functions of the
   * deal.II library
   */
  void
  project_boundary_values(const dealii::DoFHandler<dim, spacedim> &dof_handler,
                          const dealii::Quadrature<dim - 1> &      quadrature,
                          dealii::AffineConstraints<double> &constraints) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the AffineConstraints<double> .
   * It relies on the VectorTools::project_boundary_values functions of the
   * deal.II library
   */
  void
  project_boundary_values(const dealii::Mapping<dim, spacedim> &   mapping,
                          const dealii::DoFHandler<dim, spacedim> &dof_handler,
                          const dealii::Quadrature<dim - 1> &      quadrature,
                          dealii::AffineConstraints<double> &constraints) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the AffineConstraints<double> .
   * It relies on the VectorTools::project_boundary_values functions of the
   * deal.II library
   */
  void
  project_boundary_values(
    const dealii::DoFHandler<dim, spacedim> &          dof_handler,
    const dealii::Quadrature<dim - 1> &                quadrature,
    std::map<dealii::types::global_dof_index, double> &projected_bv) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the AffineConstraints<double> .
   * It relies on the VectorTools::project_boundary_values functions of the
   * deal.II library
   */
  void
  project_boundary_values(
    const dealii::Mapping<dim, spacedim> &             mapping,
    const dealii::DoFHandler<dim, spacedim> &          dof_handler,
    const dealii::Quadrature<dim - 1> &                quadrature,
    std::map<dealii::types::global_dof_index, double> &projected_bv) const;


  /**
   * This function must be called in order to apply the homogeneous Dirichlet
   * boundary conditions to the normal components of the variables specified.
   *
   * The normal component can be specified adding ".N" to the variable name
   * e.g. the normal component for the vector variable "u" will be "u.N".
   *
   * It relies on the VectorTools::compute_no_normal_flux_constraints functions
   * of the deal.II library
   */
  void
  compute_no_normal_flux_constraints(
    const dealii::DoFHandler<dim, spacedim> &dof_handler,
    dealii::AffineConstraints<double> &      constraints) const;

  /**
   * This function must be called in order to apply the homogeneous Dirichlet
   * boundary conditions to the normal components of the variables specified.
   *
   * The normal component can be specified adding ".N" to the variable name
   * e.g. the normal component for the vector variable "u" will be "u.N".
   *
   * It relies on the VectorTools::compute_no_normal_flux_constraints functions
   * of the deal.II library
   */
  void
  compute_no_normal_flux_constraints(
    const dealii::DoFHandler<dim, spacedim> &dof_handler,
    const dealii::Mapping<dim, spacedim> &   mapping,
    dealii::AffineConstraints<double> &      constraints) const;

  /**
   * This function must be called in order to apply the Dirichlet
   * boundary conditions to the normal components of the variables specified.
   *
   * The normal component can be specified adding ".N" to the variable name
   * e.g. the normal component for the vector variable "u" will be "u.N".
   *
   * It relies on the VectorTools::compute_nonzero_normal_flux_constraints
   * functions of the deal.II library.
   *
   * Note that calling this functions with an expression equal to zero is
   * equivalent to call the compute_no_normal_flux_constraints function.
   */
  void
  compute_nonzero_normal_flux_constraints(
    const dealii::DoFHandler<dim, spacedim> &dof_handler,
    dealii::AffineConstraints<double> &      constraints) const;

  /**
   * This function must be called in order to apply the Dirichlet
   * boundary conditions to the normal components of the variables specified.
   *
   * The normal component can be specified adding ".N" to the variable name
   * e.g. the normal component for the vector variable "u" will be "u.N".
   *
   * It relies on the VectorTools::compute_nonzero_normal_flux_constraints
   * functions of the deal.II library.
   *
   * Note that calling this functions with an expression equal to zero is
   * equivalent to call the compute_no_normal_flux_constraints function.
   */
  void
  compute_nonzero_normal_flux_constraints(
    const dealii::DoFHandler<dim, spacedim> &dof_handler,
    const dealii::Mapping<dim, spacedim> &   mapping,
    dealii::AffineConstraints<double> &      constraints) const;

private:
  /**
   * Number of components of the underlying Function objects.
   */
  const unsigned int n_components;
};

D2K_NAMESPACE_CLOSE

#endif
