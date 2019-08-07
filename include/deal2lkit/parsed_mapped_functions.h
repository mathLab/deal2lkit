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

#ifndef d2k_parsed_mapped_function_h
#define d2k_parsed_mapped_function_h

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parsed_function.h>

#include <deal.II/fe/component_mask.h>

#include <deal2lkit/config.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/utilities.h>

#include <algorithm>
#include <map>



D2K_NAMESPACE_OPEN

/**
 * ParsedMappedFunctions object.  It allows you to set a mapped
 * functions, i.e., parsed functions acting on specified id
 * (boundary_id, material_id,etc..) and on specified components.
 *
 * Dirichlet Boundary conditions, Neumann boundary conditions
 * and forcing terms can be easily handled with this class.
 *
 * A typical usage of this class is the following
 *
 *
 * @code
 * // the following variables are set just for the sake of completeness
 * unsigned int spacedim = 2;
 * unsigned int n_components = 3; // number of components of the problem
 *
 * // create parsed_mapped_functions object
 *
 * ParsedMappedFunctions<spacedim>
 *    parsed_mapped_functions("Forcing terms", // name for the section of the
 * Parameter Handler to use n_components, // The number of components of the
 * function "u,u,p",         // names of known components that can be used
 * instead of component numbers "0=u % 1=2 % 6=ALL", // boundary_id =
 * component;other_component % other_id = comp; other_comp "0=x;y;0 % 1=0;0;0 %
 * 6=y*k;0;k", // boundary_id = expression % other_id = other_expression "k=1");
 * // list of constants that can be used in the above epressions
 * ...
 *
 * unsigned int id = cell->material_id();
 *
 * std::vector<double> fs(n_q_points);
 *
 * parsed_mapped_functions.get_mapped_function(id)->value_list(fe_values.get_quadrature_points(),
 * fs);
 *
 * @endcode
 */

template <int spacedim>
class ParsedMappedFunctions : public deal2lkit::ParameterAcceptor
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
   * - the default names of known components that can be used instead of
   * component numbers
   *
   * - a list of ids and components where the boundary conditions must be
   * applied, which is a string with the following pattern boundary_id =
   * component;other_component % other_id = comp; other_comp
   *
   * - a list of ids and expressions defined over the ids
   * (if this string is left empty, ZeroFunction is imposed on the above
   * specified ids and components)
   *
   * - list of constants that can be used in the above epressions
   *
   */
  ParsedMappedFunctions(const std::string & name         = "Mapped Functions",
                        const unsigned int &n_components = 1,
                        const std::string & component_names       = "",
                        const std::string & default_id_components = "0=ALL",
                        const std::string & default_id_functions  = "",
                        const std::string & default_constants     = "");

  /**
   * return a shared_ptr to the ParsedFunction corresponding to the given id
   */
  shared_ptr<dealii::Functions::ParsedFunction<spacedim>>
  get_mapped_function(const unsigned int &id) const;

  /**
   * return a shared_ptr to the ParsedFunction corresponding to the given id
   * the function has spacedim components
   */
  shared_ptr<dealii::Functions::ParsedFunction<spacedim>>
  get_mapped_normal_function(const unsigned int &id,
                             const unsigned int &fcv) const;

  /**
   * return the ComponentMask corresponding to the given id
   */
  dealii::ComponentMask
  get_mapped_mask(const unsigned int &id) const;

  /**
   * return the list of the mapped ids
   */
  std::vector<unsigned int>
  get_mapped_ids() const;

  /**
   * return the list of the mapped ids for which normal components have been set
   */
  std::vector<unsigned int>
  get_mapped_normal_ids() const;


  /**
   * declare_parameters is inherithed by ParameterAcceptor
   */
  virtual void
  declare_parameters(dealii::ParameterHandler &prm);

  /**
   * parse_parameters_call_back is inherithed by ParameterAcceptor
   */
  virtual void
  parse_parameters_call_back();

  /**
   * return true if there is a function that acts on the passed id
   */
  bool
  acts_on_id(unsigned int &id) const;

  /**
   * set time equal to t for all the mapped functions
   */
  void
  set_time(const double &t);


  /// Mismatch between the number of ids set in 'IDs and component masks'
  /// and 'IDs and expressions'
  DeclException2(ExcIdsMismatch,
                 unsigned int,
                 unsigned int,
                 << "The number of ids specified in the field "
                 << "'IDs and component masks' (" << arg1 << ") "
                 << "must match the number of ids "
                 << "set in the field 'IDs and expressions' (" << arg2 << ").");

  /// An entry with this id does not exist in this object.
  DeclException1(ExcIdNotFound,
                 unsigned int,
                 << "No entry with the id " << arg1 << " exists.");

  /// No component mask are defined on this id.
  DeclException1(ExcIdNotMatch,
                 unsigned int,
                 << "No component mask associated to the id " << arg1
                 << " are defined.");

  /// Wrong number of component mask is defined on this id.
  DeclException3(ExcWrongComponent,
                 unsigned int,
                 unsigned int,
                 unsigned int,
                 << "At id " << arg1
                 << ", wrong component number has been used: " << arg2
                 << " is not in the range [0, " << arg3 << ").");
  /// Wrong variable name is defined on this id.
  DeclException3(ExcWrongVariable,
                 unsigned int,
                 std::string,
                 std::vector<std::string>,
                 << "At id " << arg1 << ", wrong variabile name has been used: "
                 << arg2 << " does not belong to the knwon variables: "
                 << print(unique(arg3)) << ".");

protected:
  void
  split_id_components(const std::string &parsed_idcomponents);
  void
  split_id_functions(const std::string &parsed_idfunctions,
                     const std::string &constants);
  void
  add_normal_components();

  void
  set_normal_functions();

  std::string                                   name;
  std::string                                   str_id_components;
  std::string                                   str_id_functions;
  std::string                                   str_component_names;
  std::string                                   str_constants;
  std::vector<std::string>                      _component_names;
  std::vector<std::string>                      _normal_components;
  std::vector<std::string>                      _all_components;
  std::vector<unsigned int>                     ids;
  std::vector<unsigned int>                     normal_ids;
  std::map<unsigned int, dealii::ComponentMask> id_components;
  std::map<unsigned int,
           shared_ptr<dealii::Functions::ParsedFunction<spacedim>>>
    id_functions;
  std::map<unsigned int,
           std::pair<dealii::ComponentMask,
                     shared_ptr<dealii::Functions::ParsedFunction<spacedim>>>>
    mapped_functions;
  std::map<std::string, std::pair<std::vector<unsigned int>, unsigned int>>
    mapped_normal_components; // name, (std::vector<boundary ids>, first
                              // component vector)
  std::vector<std::pair<unsigned int, std::string>>
                                      normal_components; // first component vector, variable_name+"N"
  std::map<unsigned int, std::string> id_str_functions;
  std::map<std::pair<unsigned int, unsigned int>,
           shared_ptr<dealii::Functions::ParsedFunction<spacedim>>>
    _normal_functions;

  const unsigned int n_components;
};



D2K_NAMESPACE_CLOSE

#endif
