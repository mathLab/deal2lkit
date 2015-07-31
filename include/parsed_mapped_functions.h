#ifndef __dealii_sak_parsed_mapped_function_h
#define __dealii_sak_parsed_mapped_function_h

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parsed_function.h>
#include "parameter_acceptor.h"
#include "parsed_function.h"
#include <deal.II/fe/component_mask.h>
#include <algorithm>
#include <map>
#include "utilities.h"

using namespace dealii;

/**
 * ParsedMappedFunctions object.
 * It allows you to set a mapped functions, i.e., parsed functions
 * acting on specified id (boundary_id, material_id,etc..) and on
 * specified components.
 *
 * Dirichlet Boundary conditions, Neumann boundary conditions
 * and forcing terms can be easily handled with this class.
 *
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
 * ParsedMappedFunctions<spacedim,n_components>
 *    parsed_mapped_functions("Forcing terms", // name for the section of the Parameter Handler to use
 *                     "u,u,p",         // names of known components that can be used instead of component numbers
 *                     "0=u % 1=2 % 6=ALL", // boundary_id = component;other_component % other_id = comp; other_comp
 *                     "0=x;y;0 % 1=0;0;0 % 6=y*k;0;k", // boundary_id = expression % other_id = other_expression
 *                     "k=1"); // list of constants that can be used in the above epressions
 * ...
 *
 * unsigned int id = cell->material_id();
 *
 * std::vector<double> fs(n_q_points);
 *
 * parsed_mapped_functions.get_mapped_function(id)->value_list(fe_values.get_quadrature_points(), fs);
 *
 * @endcode
 */

template <int spacedim, int n_components>
class ParsedMappedFunctions : public ParameterAcceptor
{
public:
  /**
   * Constructor.
   *
   * It takes:
   * - the name for the section of the Parameter Handler to use
   *
   * - the names of known components that can be used instead of component numbers
   *
   * - a list of ids and components where the boundary conditions must be applied, which is
   *   a string with the following pattern boundary_id = component;other_component % other_id = comp; other_comp
   *
   * - a list of ids and exprssions defined over the ids
   * (if this string is left empty, ZeroFunction is imposed on the above specified ids and components)
   *
   * - list of constants that can be used in the above epressions
   *
   */
  ParsedMappedFunctions  (const std::string &name = "Mapped Functions",
                          const std::string &component_names = "",
                          const std::string &default_id_components = "0=ALL",
                          const std::string &default_id_functions = "",
                          const std::string &default_constants = "");

  /**
   * return a shared_ptr to the ParsedFunction corresponding to the given id
   */
  shared_ptr<dealii::Functions::ParsedFunction<spacedim> > get_mapped_function (const unsigned int &id) const;

  /**
   * return the ComponentMask corresponding to the given id
   */
  ComponentMask get_mapped_mask (const unsigned int &id) const;

  /**
   * return the list of the mapped ids
   */
  std::vector<unsigned int> get_mapped_ids() const;

  /**
   * return the list of the mapped ids for which normal components have been set
   */
  std::vector<unsigned int> get_mapped_normal_ids() const;


  /**
   * declare_parameters is inherithed by ParameterAcceptor
   */
  virtual void declare_parameters (ParameterHandler &prm);

  /**
   * parse_parameters_call_back is inherithed by ParameterAcceptor
   */
  virtual void parse_parameters_call_back ();

  /**
   * return true if there is a function that acts on the passed id
   */
  inline bool acts_on_id (unsigned int &id) const;

  /**
   * set time equal to t for all the mapped functions
   */
  void set_time (const double &t);


  /// An entry with this id does not exist in this object.
  DeclException1(ExcIdNotFound, unsigned int,
                 << "No entry with the id " << arg1 << " exists.");

  /// No component mask are defined on this id.
  DeclException1(ExcIdNotMatch, unsigned int,
                 << "No component mask associated to the id " << arg1 << " are defined.");

  /// Wrong number of component mask is defined on this id.
  DeclException3(ExcWrongComponent, unsigned int, unsigned int, unsigned int,
                 << "At id " << arg1
                 << ", wrong component number has been used: "
                 <<  arg2 << " is not in the range [0, "
                 << arg3 <<").");
  /// Wrong variable name is defined on this id.
  DeclException3(ExcWrongVariable, unsigned int, std::string, std::vector<std::string>,
                 << "At id " << arg1
                 << ", wrong variabile name has been used: "
                 <<  arg2 << " does not belong to the knwon variables: "
                 << print(unique(arg3)) <<".");

protected:

  void split_id_components(const std::string &parsed_idcomponents);
  void split_id_functions(const std::string &parsed_idfunctions,
                          const std::string &constants);
  void add_normal_components();


  std::string name;
  std::string str_id_components;
  std::string str_id_functions;
  std::string str_component_names;
  std::string str_constants;
  std::vector<std::string> _component_names;
  std::vector<std::string> _normal_components;
  std::vector<std::string> _all_components;
  std::vector<unsigned int> ids;
  std::vector<unsigned int> normal_ids;
  std::map<unsigned int, ComponentMask> id_components;
  std::map<unsigned int, shared_ptr<dealii::Functions::ParsedFunction<spacedim> > > id_functions;
  std::map<unsigned int, shared_ptr<dealii::Functions::ParsedFunction<spacedim> > > normal_id_functions;
  std::map<unsigned int, std::pair<ComponentMask, shared_ptr<dealii::Functions::ParsedFunction<spacedim> > > > mapped_functions;
  std::map<std::string, std::pair<std::vector<unsigned int>, unsigned int > > mapped_normal_components; // name, (std::vector<boundary ids>, firs component vector)
  std::vector<std::pair<unsigned int, std::string> > normal_components; // first component vector, variable_name+"N"

};


#endif
