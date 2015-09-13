#ifndef __dealii_sak_parsed_dirichlet_bcs_h
#define __dealii_sak_parsed_dirichlet_bcs_h

#include <deal2lkit/parameter_acceptor.h>
#include <deal.II/base/exceptions.h>
#include <deal2lkit/parsed_mapped_functions.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;
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
 * ParsedDirichletBCs<dim,spacedim,n_components>
 *    parsed_dirichlet("Dirichlet BCs", // name for the section of the Parameter Handler to use
 *                     "u,u,p",         // names of known components that can be used instead of component numbers
 *                     "0=u % 1=2 % 3=u.N % 6=ALL", // boundary_id = component;other_component % other_id = comp; other_comp
 *                     "0=x;y;0 % 1=0;0;0 % 3=x;2;0 % 6=y*k;0;k", // boundary_id = expression % other_id = other_expression
 *                     "k=1"); // list of constants that can be used in the above epressions
 * ...
 * QGauss<dim-1> quadrature(2);
 * MappingQ1<dim> mapping;
 * std::map<types::global_dof_index,double> dirichlet_dofs;
 * ConstraintMatrix constraints;
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

template <int dim, int spacedim, int n_components>
class ParsedDirichletBCs : public ParsedMappedFunctions<spacedim,n_components>
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
   * (if this string is left empty, homogeneous boundary conditions are imposed on the above specified ids and components)
   *
   * - list of constants that can be used in the above epressions
   *
   */
  ParsedDirichletBCs (const std::string &name = "Dirichlet BCs",
                      const std::string &component_names = "",
                      const std::string &default_id_components = "0=ALL",
                      const std::string &default_id_functions = "",
                      const std::string &default_constants = "");

  /**
   * these method calls the method of the Parent class
   */
  virtual void declare_parameters (ParameterHandler &prm);

  /**
   * these method calls the method of the Parent class
   */
  virtual void parse_parameters_call_back ();

  /**
   * This function must be called in order to apply the boundary conditions
   * to the ConstraintMatrix.
   * It relies on the VectorTools::interpolate_boundary_values functions of the
   * deal.II library
   */
  void interpolate_boundary_values (const DoFHandler<dim,spacedim> &dof_handler,
                                    ConstraintMatrix &constraints) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the ConstraintMatrix.
   * It relies on the VectorTools::interpolate_boundary_values functions of the
   * deal.II library
   */
  void interpolate_boundary_values (const Mapping<dim,spacedim> &mapping,
                                    const DoFHandler<dim,spacedim> &dof_handler,
                                    ConstraintMatrix &constraints) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the ConstraintMatrix.
   * It relies on the VectorTools::interpolate_boundary_values functions of the
   * deal.II library
   */
  void interpolate_boundary_values (const DoFHandler<dim,spacedim> &dof_handler,
                                    std::map<types::global_dof_index,double> &d_dofs) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the ConstraintMatrix.
   * It relies on the VectorTools::interpolate_boundary_values functions of the
   * deal.II library
   */
  void interpolate_boundary_values (const Mapping<dim,spacedim> &mapping,
                                    const DoFHandler<dim,spacedim> &dof_handler,
                                    std::map<types::global_dof_index,double> &d_dofs) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the ConstraintMatrix.
   * It relies on the VectorTools::project_boundary_values functions of the
   * deal.II library
   */
  void project_boundary_values (const DoFHandler<dim,spacedim> &dof_handler,
                                const Quadrature<dim-1> &quadrature,
                                ConstraintMatrix &constraints) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the ConstraintMatrix.
   * It relies on the VectorTools::project_boundary_values functions of the
   * deal.II library
   */
  void project_boundary_values (const Mapping<dim,spacedim> &mapping,
                                const DoFHandler<dim,spacedim> &dof_handler,
                                const Quadrature<dim-1> &quadrature,
                                ConstraintMatrix &constraints) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the ConstraintMatrix.
   * It relies on the VectorTools::project_boundary_values functions of the
   * deal.II library
   */
  void project_boundary_values (const DoFHandler<dim,spacedim> &dof_handler,
                                const Quadrature<dim-1> &quadrature,
                                std::map<types::global_dof_index,double> &projected_bv) const;

  /**
   * This function must be called in order to apply the boundary conditions
   * to the ConstraintMatrix.
   * It relies on the VectorTools::project_boundary_values functions of the
   * deal.II library
   */
  void project_boundary_values (const Mapping<dim,spacedim> &mapping,
                                const DoFHandler<dim,spacedim> &dof_handler,
                                const Quadrature<dim-1> &quadrature,
                                std::map<types::global_dof_index,double> &projected_bv) const;


  /**
   * This function must be called in order to apply the homogeneous Dirichlet
   * boundary conditions to the normal components of the variables specified.
   *
   * The normal component can be specified adding ".N" to the variable name
   * e.g. the normal component for the vector variable "u" will be "u.N".
   *
   * It relies on the VectorTools::compute_no_normal_flux_constraints functions of the
   * deal.II library
   */
  void compute_no_normal_flux_constraints (const DoFHandler<dim,spacedim> &dof_handler,
                                           ConstraintMatrix &constraints) const;
  /**
   * This function must be called in order to apply the homogeneous Dirichlet
   * boundary conditions to the normal components of the variables specified.
   *
   * The normal component can be specified adding ".N" to the variable name
   * e.g. the normal component for the vector variable "u" will be "u.N".
   *
   * It relies on the VectorTools::compute_no_normal_flux_constraints functions of the
   * deal.II library
   */
  void compute_no_normal_flux_constraints (const DoFHandler<dim,spacedim> &dof_handler,
                                           const Mapping< dim, spacedim > &mapping,
                                           ConstraintMatrix &constraints) const;
  /**
   * This function must be called in order to apply the Dirichlet
   * boundary conditions to the normal components of the variables specified.
   *
   * The normal component can be specified adding ".N" to the variable name
   * e.g. the normal component for the vector variable "u" will be "u.N".
   *
   * It relies on the VectorTools::compute_nonzero_normal_flux_constraints functions of the
   * deal.II library.
   *
   * Note that calling this functions with an expression equal to zero is equivalent to call
   * the compute_no_normal_flux_constraints function.
   */
  void compute_nonzero_normal_flux_constraints (const DoFHandler<dim,spacedim> &dof_handler,
                                                ConstraintMatrix &constraints) const;

  /**
   * This function must be called in order to apply the Dirichlet
   * boundary conditions to the normal components of the variables specified.
   *
   * The normal component can be specified adding ".N" to the variable name
   * e.g. the normal component for the vector variable "u" will be "u.N".
   *
   * It relies on the VectorTools::compute_nonzero_normal_flux_constraints functions of the
   * deal.II library.
   *
   * Note that calling this functions with an expression equal to zero is equivalent to call
   * the compute_no_normal_flux_constraints function.
   */
  void compute_nonzero_normal_flux_constraints (const DoFHandler<dim,spacedim> &dof_handler,
                                                const Mapping< dim, spacedim > &mapping,
                                                ConstraintMatrix &constraints) const;

};


template <int dim, int spacedim, int n_components>
ParsedDirichletBCs<dim,spacedim,n_components>::
ParsedDirichletBCs (const std::string &parsed_name,
                    const std::string &parsed_component_names,
                    const std::string &parsed_id_components,
                    const std::string &parsed_id_functions,
                    const std::string &parsed_constants)
  :

  ParsedMappedFunctions<spacedim,n_components>(parsed_name,
                                               parsed_component_names,
                                               parsed_id_components,
                                               parsed_id_functions,
                                               parsed_constants)
{}

template<int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::declare_parameters (ParameterHandler &prm)
{
  ParsedMappedFunctions<spacedim,n_components>::declare_parameters(prm);
}

template<int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::parse_parameters_call_back ()
{
  ParsedMappedFunctions<spacedim,n_components>::parse_parameters_call_back();
}

template<int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::interpolate_boundary_values(const DoFHandler<dim,spacedim> &dof_handler,
    ConstraintMatrix &constraints) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i=0; i<ids.size(); ++i)
    VectorTools::interpolate_boundary_values (dof_handler,
                                              ids[i],
                                              *(this->get_mapped_function(ids[i])),
                                              constraints,
                                              this->get_mapped_mask(ids[i]));
}

template<int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::interpolate_boundary_values(const Mapping<dim,spacedim> &mapping,
    const DoFHandler<dim,spacedim> &dof_handler,
    ConstraintMatrix &constraints) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i=0; i<ids.size(); ++i)
    VectorTools::interpolate_boundary_values (mapping,
                                              dof_handler,
                                              ids[i],
                                              *(this->get_mapped_function(ids[i])),
                                              constraints,
                                              this->get_mapped_mask(ids[i]));
}

template<int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::interpolate_boundary_values( const DoFHandler<dim,spacedim> &dof_handler,
    std::map<types::global_dof_index,double> &d_dofs) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i=0; i<ids.size(); ++i)
    VectorTools::interpolate_boundary_values (dof_handler,
                                              ids[i],
                                              *(this->get_mapped_function(ids[i])),
                                              d_dofs,
                                              this->get_mapped_mask(ids[i]));
}

template<int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::interpolate_boundary_values(const Mapping<dim,spacedim> &mapping,
    const DoFHandler<dim,spacedim> &dof_handler,
    std::map<types::global_dof_index,double> &d_dofs) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i=0; i<ids.size(); ++i)
    VectorTools::interpolate_boundary_values (mapping,
                                              dof_handler,
                                              ids[i],
                                              *(this->get_mapped_function(ids[i])),
                                              d_dofs,
                                              this->get_mapped_mask(ids[i]));
}

template <int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::project_boundary_values (const Mapping<dim,spacedim> &mapping,
    const DoFHandler<dim,spacedim> &dof_handler,
    const Quadrature<dim-1> &quadrature,
    ConstraintMatrix &constraints) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i=0; i<ids.size(); ++i)
    {
      typename FunctionMap<spacedim>::type boundary_map;

      Function<spacedim> *f;
      f = &(*(this->get_mapped_function(ids[i])));
      boundary_map[ids[i]] = f;

      // from bool to int
      std::vector<unsigned int> component_mapping;
      if (spacedim > 1) // component_mapping is not supported in deal.II for dim==1
        for (unsigned int j=0; j<component_mapping.size(); ++j)
          component_mapping.push_back(this->get_mapped_mask(ids[i])[j]);

      VectorTools::project_boundary_values(mapping,
                                           dof_handler,
                                           boundary_map,
                                           quadrature,
                                           constraints,
                                           component_mapping);
    }
}

template <int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::project_boundary_values (const DoFHandler<dim,spacedim> &dof_handler,
    const Quadrature<dim-1> &quadrature,
    ConstraintMatrix &constraints) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i=0; i<ids.size(); ++i)
    {
      typename FunctionMap<spacedim>::type boundary_map;

      Function<spacedim> *f;
      f = &(*(this->get_mapped_function(ids[i])));
      boundary_map[ids[i]] = f;

      // from bool to int
      std::vector<unsigned int> component_mapping;
      if (spacedim > 1) // component_mapping is not supported in deal.II for dim==1
        for (unsigned int j=0; j<component_mapping.size(); ++j)
          component_mapping.push_back(this->get_mapped_mask(ids[i])[j]);

      VectorTools::project_boundary_values(dof_handler,
                                           boundary_map,
                                           quadrature,
                                           constraints,
                                           component_mapping);
    }
}

template <int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::project_boundary_values (const Mapping<dim,spacedim> &mapping,
    const DoFHandler<dim,spacedim> &dof_handler,
    const Quadrature<dim-1> &quadrature,
    std::map<types::global_dof_index,double> &projected_bv) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i=0; i<ids.size(); ++i)
    {
      typename FunctionMap<spacedim>::type boundary_map;

      Function<spacedim> *f;
      f = &(*(this->get_mapped_function(ids[i])));
      boundary_map[ids[i]] = f;

      // from bool to int
      std::vector<unsigned int> component_mapping;
      if (spacedim > 1) // component_mapping is not supported in deal.II for dim==1
        for (unsigned int j=0; j<n_components; ++j)
          component_mapping.push_back(this->get_mapped_mask(ids[i])[j]);

      VectorTools::project_boundary_values(mapping,
                                           dof_handler,
                                           boundary_map,
                                           quadrature,
                                           projected_bv,
                                           component_mapping);

    }
}


template <int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::project_boundary_values (const DoFHandler<dim,spacedim> &dof_handler,
    const Quadrature<dim-1> &quadrature,
    std::map<types::global_dof_index,double> &projected_bv) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i=0; i<ids.size(); ++i)
    {
      typename FunctionMap<spacedim>::type boundary_map;

      Function<spacedim> *f;
      f = &(*(this->get_mapped_function(ids[i])));
      boundary_map[ids[i]] = f;

      // from bool to int
      std::vector<unsigned int> component_mapping;
      if (spacedim > 1) // component_mapping is not supported in deal.II for dim==1
        for (unsigned int j=0; j<n_components; ++j)
          component_mapping.push_back(this->get_mapped_mask(ids[i])[j]);

      VectorTools::project_boundary_values(dof_handler,
                                           boundary_map,
                                           quadrature,
                                           projected_bv,
                                           component_mapping);
    }
}

template <int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::compute_no_normal_flux_constraints(const DoFHandler<dim,spacedim> &dof_handler,
    ConstraintMatrix &constraints) const
{
  std::set<types::boundary_id> no_normal_flux_boundaries;

  typedef std::map<std::string, std::pair<std::vector<unsigned int>, unsigned int > >::const_iterator it_type;

  for (it_type it=this->mapped_normal_components.begin(); it != this->mapped_normal_components.end(); ++it)
    {
      std::vector<unsigned int> normal_ids = (it->second).first;

      for (unsigned int i=0; i<normal_ids.size(); ++i)
        no_normal_flux_boundaries.insert(normal_ids[i]);

      VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                      (it->second).second, // unsigned int first component vector
                                                      no_normal_flux_boundaries,
                                                      constraints);
    }
}

template <int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::compute_no_normal_flux_constraints(const DoFHandler<dim,spacedim> &dof_handler,
    const Mapping< dim, spacedim > &mapping,
    ConstraintMatrix &constraints) const
{
  std::set<types::boundary_id> no_normal_flux_boundaries;

  typedef std::map<std::string, std::pair<std::vector<unsigned int>, unsigned int > >::const_iterator it_type;

  for (it_type it=this->mapped_normal_components.begin(); it != this->mapped_normal_components.end(); ++it)
    {
      std::vector<unsigned int> normal_ids = (it->second).first;

      for (unsigned int i=0; i<normal_ids.size(); ++i)
        no_normal_flux_boundaries.insert(normal_ids[i]);

      VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                      (it->second).second, // unsigned int first component vector
                                                      no_normal_flux_boundaries,
                                                      constraints,
                                                      mapping);
    }
}


template <int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::compute_nonzero_normal_flux_constraints(const DoFHandler<dim,spacedim> &dof_handler,
    ConstraintMatrix &constraints) const
{
  std::set<types::boundary_id> no_normal_flux_boundaries;

  typedef std::map<std::string, std::pair<std::vector<unsigned int>, unsigned int > >::const_iterator it_type;

  for (it_type it=this->mapped_normal_components.begin(); it != this->mapped_normal_components.end(); ++it)
    {
      typename FunctionMap<spacedim>::type boundary_map;

      std::vector<unsigned int> normal_ids = (it->second).first;
      unsigned int fcv = (it->second).second; // unsigned int first component vector

      for (unsigned int i=0; i<normal_ids.size(); ++i)
        {
          Function<spacedim> *f;
          f = &(*(this->get_mapped_normal_function(normal_ids[i], fcv)));
          boundary_map[normal_ids[i]] = f;
          no_normal_flux_boundaries.insert(normal_ids[i]);
        }

      VectorTools::compute_nonzero_normal_flux_constraints(dof_handler,
                                                           fcv,
                                                           no_normal_flux_boundaries,
                                                           boundary_map,
                                                           constraints);
    }
}

template <int dim, int spacedim, int n_components>
void ParsedDirichletBCs<dim,spacedim,n_components>::compute_nonzero_normal_flux_constraints(const DoFHandler<dim,spacedim> &dof_handler,
    const Mapping< dim, spacedim > &mapping,
    ConstraintMatrix &constraints) const
{
  std::set<types::boundary_id> no_normal_flux_boundaries;

  typedef std::map<std::string, std::pair<std::vector<unsigned int>, unsigned int > >::const_iterator it_type;

  for (it_type it=this->mapped_normal_components.begin(); it != this->mapped_normal_components.end(); ++it)
    {
      typename FunctionMap<spacedim>::type boundary_map;

      std::vector<unsigned int> normal_ids = (it->second).first;
      unsigned int fcv = (it->second).second; // unsigned int first component vector

      for (unsigned int i=0; i<normal_ids.size(); ++i)
        {
          Function<spacedim> *f;

          f = &(*(this->get_mapped_normal_function(normal_ids[i], fcv)));
          boundary_map[normal_ids[i]] = f;
          no_normal_flux_boundaries.insert(normal_ids[i]);
        }

      VectorTools::compute_nonzero_normal_flux_constraints(dof_handler,
                                                           fcv,
                                                           no_normal_flux_boundaries,
                                                           boundary_map,
                                                           constraints,
                                                           mapping);
    }
}


#endif
