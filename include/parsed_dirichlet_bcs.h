#ifndef __dealii_sak_parsed_dirichlet_bcs_h
#define __dealii_sak_parsed_dirichlet_bcs_h

#include "parameter_acceptor.h"
#include <deal.II/base/exceptions.h>
#include "parsed_mapped_functions.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

template <int dim, int spacedim, int n_components>
class ParsedDirichletBCs : public ParsedMappedFunctions<spacedim,n_components>
{
public:
  ParsedDirichletBCs (const std::string &name = "Dirichlet BCs",
                      const std::string &component_names = "",
                      const std::string &default_id_components = "0=ALL",
                      const std::string &default_id_functions = "",
                      const std::string &default_constants = "");
  virtual void declare_parameters (ParameterHandler &prm);
  virtual void parse_parameters_call_back ();

  void interpolate_boundary_values (const DoFHandler<dim,spacedim> &dof_handler,
                                    ConstraintMatrix &constraints) const;

  void interpolate_boundary_values (const Mapping<dim,spacedim> &mapping,
                                    const DoFHandler<dim,spacedim> &dof_handler,
                                    ConstraintMatrix &constraints) const;

  void interpolate_boundary_values (const DoFHandler<dim,spacedim> &dof_handler,
                                    std::map<types::global_dof_index,double> &d_dofs) const;
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


#endif
