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

#include <deal2lkit/parsed_dirichlet_bcs.h>

using namespace dealii;

D2K_NAMESPACE_OPEN

template <int dim, int spacedim>
ParsedDirichletBCs<dim, spacedim>::ParsedDirichletBCs(
  const std::string & parsed_name,
  const unsigned int &n_components,
  const std::string & parsed_component_names,
  const std::string & parsed_id_components,
  const std::string & parsed_id_functions,
  const std::string & parsed_constants)
  :

  ParsedMappedFunctions<spacedim>(parsed_name,
                                  n_components,
                                  parsed_component_names,
                                  parsed_id_components,
                                  parsed_id_functions,
                                  parsed_constants)
  , n_components(n_components)
{}

template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::declare_parameters(ParameterHandler &prm)
{
  ParsedMappedFunctions<spacedim>::declare_parameters(prm);
}

template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::parse_parameters_call_back()
{
  ParsedMappedFunctions<spacedim>::parse_parameters_call_back();
}

template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::interpolate_boundary_values(
  const DoFHandler<dim, spacedim> &dof_handler,
  ConstraintMatrix &               constraints) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i = 0; i < ids.size(); ++i)
    VectorTools::interpolate_boundary_values(dof_handler,
                                             ids[i],
                                             *(this->get_mapped_function(
                                               ids[i])),
                                             constraints,
                                             this->get_mapped_mask(ids[i]));
}

template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::interpolate_boundary_values(
  const Mapping<dim, spacedim> &   mapping,
  const DoFHandler<dim, spacedim> &dof_handler,
  ConstraintMatrix &               constraints) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i = 0; i < ids.size(); ++i)
    VectorTools::interpolate_boundary_values(mapping,
                                             dof_handler,
                                             ids[i],
                                             *(this->get_mapped_function(
                                               ids[i])),
                                             constraints,
                                             this->get_mapped_mask(ids[i]));
}

template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::interpolate_boundary_values(
  const DoFHandler<dim, spacedim> &          dof_handler,
  std::map<types::global_dof_index, double> &d_dofs) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i = 0; i < ids.size(); ++i)
    VectorTools::interpolate_boundary_values(dof_handler,
                                             ids[i],
                                             *(this->get_mapped_function(
                                               ids[i])),
                                             d_dofs,
                                             this->get_mapped_mask(ids[i]));
}

template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::interpolate_boundary_values(
  const Mapping<dim, spacedim> &             mapping,
  const DoFHandler<dim, spacedim> &          dof_handler,
  std::map<types::global_dof_index, double> &d_dofs) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i = 0; i < ids.size(); ++i)
    VectorTools::interpolate_boundary_values(mapping,
                                             dof_handler,
                                             ids[i],
                                             *(this->get_mapped_function(
                                               ids[i])),
                                             d_dofs,
                                             this->get_mapped_mask(ids[i]));
}


// [TODO] Fix this in deal.II.

template <>
void ParsedDirichletBCs<1, 2>::project_boundary_values(DoFHandler<1, 2> const &,
                                                       Quadrature<0> const &,
                                                       ConstraintMatrix &) const
{
  ExcImpossibleInDim(1);
}

template <>
void ParsedDirichletBCs<1, 3>::project_boundary_values(DoFHandler<1, 3> const &,
                                                       Quadrature<0> const &,
                                                       ConstraintMatrix &) const
{
  ExcImpossibleInDim(1);
}

template <>
void ParsedDirichletBCs<2, 3>::project_boundary_values(DoFHandler<2, 3> const &,
                                                       Quadrature<1> const &,
                                                       ConstraintMatrix &) const
{
  ExcNotImplemented();
}


template <>
void
ParsedDirichletBCs<1, 2>::project_boundary_values(const Mapping<1, 2> &,
                                                  DoFHandler<1, 2> const &,
                                                  Quadrature<0> const &,
                                                  ConstraintMatrix &) const
{
  ExcImpossibleInDim(1);
}

template <>
void
ParsedDirichletBCs<1, 3>::project_boundary_values(const Mapping<1, 3> &,
                                                  DoFHandler<1, 3> const &,
                                                  Quadrature<0> const &,
                                                  ConstraintMatrix &) const
{
  ExcImpossibleInDim(1);
}



template <>
void
ParsedDirichletBCs<2, 3>::project_boundary_values(const Mapping<2, 3> &,
                                                  DoFHandler<2, 3> const &,
                                                  Quadrature<1> const &,
                                                  ConstraintMatrix &) const
{
  ExcNotImplemented();
}



template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::project_boundary_values(
  const Mapping<dim, spacedim> &   mapping,
  const DoFHandler<dim, spacedim> &dof_handler,
  const Quadrature<dim - 1> &      quadrature,
  ConstraintMatrix &               constraints) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i = 0; i < ids.size(); ++i)
    {
      typename FunctionMap<spacedim>::type boundary_map;

      Function<spacedim> *f;
      f                    = &(*(this->get_mapped_function(ids[i])));
      boundary_map[ids[i]] = f;

      // from bool to int
      std::vector<unsigned int> component_mapping;
      if (spacedim >
          1) // component_mapping is not supported in deal.II for dim==1
        for (unsigned int j = 0; j < component_mapping.size(); ++j)
          component_mapping.push_back(this->get_mapped_mask(ids[i])[j]);

      VectorTools::project_boundary_values(mapping,
                                           dof_handler,
                                           boundary_map,
                                           quadrature,
                                           constraints,
                                           component_mapping);
    }
}

template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::project_boundary_values(
  const DoFHandler<dim, spacedim> &dof_handler,
  const Quadrature<dim - 1> &      quadrature,
  ConstraintMatrix &               constraints) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i = 0; i < ids.size(); ++i)
    {
      typename FunctionMap<spacedim>::type boundary_map;

      Function<spacedim> *f;
      f                    = &(*(this->get_mapped_function(ids[i])));
      boundary_map[ids[i]] = f;

      // from bool to int
      std::vector<unsigned int> component_mapping;
      if (spacedim >
          1) // component_mapping is not supported in deal.II for dim==1
        for (unsigned int j = 0; j < component_mapping.size(); ++j)
          component_mapping.push_back(this->get_mapped_mask(ids[i])[j]);

      VectorTools::project_boundary_values(
        dof_handler, boundary_map, quadrature, constraints, component_mapping);
    }
}

template <>
void ParsedDirichletBCs<1, 2>::project_boundary_values(
  Mapping<1, 2> const &,
  DoFHandler<1, 2> const &,
  Quadrature<0> const &,
  std::map<types::global_dof_index, double> &) const
{
  Assert(false, ExcImpossibleInDim(1));
}

template <>
void ParsedDirichletBCs<1, 3>::project_boundary_values(
  Mapping<1, 3> const &,
  DoFHandler<1, 3> const &,
  Quadrature<0> const &,
  std::map<types::global_dof_index, double> &) const
{
  Assert(false, ExcImpossibleInDim(1));
}



template <>
void ParsedDirichletBCs<2, 3>::project_boundary_values(
  Mapping<2, 3> const &,
  DoFHandler<2, 3> const &,
  Quadrature<1> const &,
  std::map<types::global_dof_index, double> &) const
{
  Assert(false, ExcImpossibleInDim(2));
}


template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::project_boundary_values(
  const Mapping<dim, spacedim> &             mapping,
  const DoFHandler<dim, spacedim> &          dof_handler,
  const Quadrature<dim - 1> &                quadrature,
  std::map<types::global_dof_index, double> &projected_bv) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i = 0; i < ids.size(); ++i)
    {
      typename FunctionMap<spacedim>::type boundary_map;

      Function<spacedim> *f;
      f                    = &(*(this->get_mapped_function(ids[i])));
      boundary_map[ids[i]] = f;

      // from bool to int
      std::vector<unsigned int> component_mapping;
      if (spacedim >
          1) // component_mapping is not supported in deal.II for dim==1
        for (unsigned int j = 0; j < n_components; ++j)
          component_mapping.push_back(this->get_mapped_mask(ids[i])[j]);

      VectorTools::project_boundary_values(mapping,
                                           dof_handler,
                                           boundary_map,
                                           quadrature,
                                           projected_bv,
                                           component_mapping);
    }
}


template <>
void
ParsedDirichletBCs<1, 2>::project_boundary_values(
  const DoFHandler<1, 2> &,
  const Quadrature<0> &,
  std::map<types::global_dof_index, double> &) const
{
  Assert(false, ExcImpossibleInDim(1));
}


template <>
void
ParsedDirichletBCs<1, 3>::project_boundary_values(
  const DoFHandler<1, 3> &,
  const Quadrature<0> &,
  std::map<types::global_dof_index, double> &) const
{
  Assert(false, ExcImpossibleInDim(1));
}


template <>
void
ParsedDirichletBCs<2, 3>::project_boundary_values(
  const DoFHandler<2, 3> &,
  const Quadrature<1> &,
  std::map<types::global_dof_index, double> &) const
{
  Assert(false, ExcNotImplemented());
}


template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::project_boundary_values(
  const DoFHandler<dim, spacedim> &          dof_handler,
  const Quadrature<dim - 1> &                quadrature,
  std::map<types::global_dof_index, double> &projected_bv) const
{
  std::vector<unsigned int> ids = this->get_mapped_ids();
  for (unsigned int i = 0; i < ids.size(); ++i)
    {
      typename FunctionMap<spacedim>::type boundary_map;

      Function<spacedim> *f;
      f                    = &(*(this->get_mapped_function(ids[i])));
      boundary_map[ids[i]] = f;

      // from bool to int
      std::vector<unsigned int> component_mapping;
      if (spacedim >
          1) // component_mapping is not supported in deal.II for dim==1
        for (unsigned int j = 0; j < n_components; ++j)
          component_mapping.push_back(this->get_mapped_mask(ids[i])[j]);

      VectorTools::project_boundary_values(
        dof_handler, boundary_map, quadrature, projected_bv, component_mapping);
    }
}



template <>
void
ParsedDirichletBCs<1, 1>::compute_no_normal_flux_constraints(
  const DoFHandler<1, 1> &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(1));
}

template <>
void
ParsedDirichletBCs<1, 1>::compute_no_normal_flux_constraints(
  const DoFHandler<1, 1> &,
  const Mapping<1, 1> &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(1));
}



template <>
void ParsedDirichletBCs<1, 1>::compute_nonzero_normal_flux_constraints(
  DoFHandler<1, 1> const &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(1));
}

template <>
void
ParsedDirichletBCs<1, 1>::compute_nonzero_normal_flux_constraints(
  const DoFHandler<1, 1> &,
  const Mapping<1, 1> &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(1));
}



template <>
void
ParsedDirichletBCs<1, 2>::compute_no_normal_flux_constraints(
  const DoFHandler<1, 2> &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(1));
}

template <>
void
ParsedDirichletBCs<1, 2>::compute_no_normal_flux_constraints(
  const DoFHandler<1, 2> &,
  const Mapping<1, 2> &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(1));
}



template <>
void ParsedDirichletBCs<1, 2>::compute_nonzero_normal_flux_constraints(
  DoFHandler<1, 2> const &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(1));
}

template <>
void
ParsedDirichletBCs<1, 2>::compute_nonzero_normal_flux_constraints(
  const DoFHandler<1, 2> &,
  const Mapping<1, 2> &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(1));
}


template <>
void
ParsedDirichletBCs<1, 3>::compute_no_normal_flux_constraints(
  const DoFHandler<1, 3> &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(1));
}

template <>
void
ParsedDirichletBCs<1, 3>::compute_no_normal_flux_constraints(
  const DoFHandler<1, 3> &,
  const Mapping<1, 3> &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(1));
}



template <>
void ParsedDirichletBCs<1, 3>::compute_nonzero_normal_flux_constraints(
  DoFHandler<1, 3> const &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(1));
}

template <>
void
ParsedDirichletBCs<1, 3>::compute_nonzero_normal_flux_constraints(
  const DoFHandler<1, 3> &,
  const Mapping<1, 3> &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(1));
}



template <>
void
ParsedDirichletBCs<2, 3>::compute_no_normal_flux_constraints(
  const DoFHandler<2, 3> &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(2));
}

template <>
void
ParsedDirichletBCs<2, 3>::compute_no_normal_flux_constraints(
  const DoFHandler<2, 3> &,
  const Mapping<2, 3> &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(2));
}



template <>
void ParsedDirichletBCs<2, 3>::compute_nonzero_normal_flux_constraints(
  DoFHandler<2, 3> const &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(2));
}

template <>
void
ParsedDirichletBCs<2, 3>::compute_nonzero_normal_flux_constraints(
  const DoFHandler<2, 3> &,
  const Mapping<2, 3> &,
  ConstraintMatrix &) const
{
  Assert(false, ExcImpossibleInDim(2));
}


template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::compute_no_normal_flux_constraints(
  const DoFHandler<dim, spacedim> &dof_handler,
  ConstraintMatrix &               constraints) const
{
  std::set<types::boundary_id> no_normal_flux_boundaries;

  typedef std::map<
    std::string,
    std::pair<std::vector<unsigned int>, unsigned int>>::const_iterator it_type;

  for (it_type it = this->mapped_normal_components.begin();
       it != this->mapped_normal_components.end();
       ++it)
    {
      std::vector<unsigned int> normal_ids = (it->second).first;

      for (unsigned int i = 0; i < normal_ids.size(); ++i)
        no_normal_flux_boundaries.insert(normal_ids[i]);

      VectorTools::compute_no_normal_flux_constraints(
        dof_handler,
        (it->second).second, // unsigned int first component vector
        no_normal_flux_boundaries,
        constraints);
    }
}

template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::compute_no_normal_flux_constraints(
  const DoFHandler<dim, spacedim> &dof_handler,
  const Mapping<dim, spacedim> &   mapping,
  ConstraintMatrix &               constraints) const
{
  std::set<types::boundary_id> no_normal_flux_boundaries;

  typedef std::map<
    std::string,
    std::pair<std::vector<unsigned int>, unsigned int>>::const_iterator it_type;

  for (it_type it = this->mapped_normal_components.begin();
       it != this->mapped_normal_components.end();
       ++it)
    {
      std::vector<unsigned int> normal_ids = (it->second).first;

      for (unsigned int i = 0; i < normal_ids.size(); ++i)
        no_normal_flux_boundaries.insert(normal_ids[i]);

      VectorTools::compute_no_normal_flux_constraints(
        dof_handler,
        (it->second).second, // unsigned int first component vector
        no_normal_flux_boundaries,
        constraints,
        mapping);
    }
}


template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::compute_nonzero_normal_flux_constraints(
  const DoFHandler<dim, spacedim> &dof_handler,
  ConstraintMatrix &               constraints) const
{
  std::set<types::boundary_id> no_normal_flux_boundaries;

  typedef std::map<
    std::string,
    std::pair<std::vector<unsigned int>, unsigned int>>::const_iterator it_type;

  for (it_type it = this->mapped_normal_components.begin();
       it != this->mapped_normal_components.end();
       ++it)
    {
      typename FunctionMap<spacedim>::type boundary_map;

      std::vector<unsigned int> normal_ids = (it->second).first;
      unsigned int              fcv =
        (it->second).second; // unsigned int first component vector

      for (unsigned int i = 0; i < normal_ids.size(); ++i)
        {
          Function<spacedim> *f;
          f = &(*(this->get_mapped_normal_function(normal_ids[i], fcv)));
          boundary_map[normal_ids[i]] = f;
          no_normal_flux_boundaries.insert(normal_ids[i]);
        }

      VectorTools::compute_nonzero_normal_flux_constraints(
        dof_handler, fcv, no_normal_flux_boundaries, boundary_map, constraints);
    }
}

template <int dim, int spacedim>
void
ParsedDirichletBCs<dim, spacedim>::compute_nonzero_normal_flux_constraints(
  const DoFHandler<dim, spacedim> &dof_handler,
  const Mapping<dim, spacedim> &   mapping,
  ConstraintMatrix &               constraints) const
{
  std::set<types::boundary_id> no_normal_flux_boundaries;

  typedef std::map<
    std::string,
    std::pair<std::vector<unsigned int>, unsigned int>>::const_iterator it_type;

  for (it_type it = this->mapped_normal_components.begin();
       it != this->mapped_normal_components.end();
       ++it)
    {
      typename FunctionMap<spacedim>::type boundary_map;

      std::vector<unsigned int> normal_ids = (it->second).first;
      unsigned int              fcv =
        (it->second).second; // unsigned int first component vector

      for (unsigned int i = 0; i < normal_ids.size(); ++i)
        {
          Function<spacedim> *f;

          f = &(*(this->get_mapped_normal_function(normal_ids[i], fcv)));
          boundary_map[normal_ids[i]] = f;
          no_normal_flux_boundaries.insert(normal_ids[i]);
        }

      VectorTools::compute_nonzero_normal_flux_constraints(
        dof_handler,
        fcv,
        no_normal_flux_boundaries,
        boundary_map,
        constraints,
        mapping);
    }
}

template class ParsedDirichletBCs<1, 1>;
template class ParsedDirichletBCs<1, 2>;
template class ParsedDirichletBCs<1, 3>;
template class ParsedDirichletBCs<2, 2>;
template class ParsedDirichletBCs<2, 3>;
template class ParsedDirichletBCs<3, 3>;

D2K_NAMESPACE_CLOSE
