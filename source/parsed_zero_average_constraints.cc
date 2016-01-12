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

#include <deal2lkit/parsed_zero_average_constraints.h>
#include <deal.II/dofs/dof_tools.h>



D2K_NAMESPACE_OPEN

template <int dim, int spacedim>
ParsedZeroAverageConstraints<dim,spacedim>::
ParsedZeroAverageConstraints(const std::string &parsed_name,
                             const unsigned int &n_components,
                             const std::string &parsed_component_names,
                             const std::string &parsed_components,
                             const std::string &parsed_boundary_components):
  ParameterAcceptor(parsed_name),
  name (parsed_name),
  str_components (parsed_components),
  str_boundary_components (parsed_boundary_components),
  str_component_names (parsed_component_names),
  mask(n_components, false),
  boundary_mask(n_components, false),
  n_components(n_components)
{}

template <int dim, int spacedim>
void ParsedZeroAverageConstraints<dim,spacedim>::parse_parameters_call_back()
{
  for (unsigned int c=0; c<components.size(); ++c)
    {
      if ((std::find(_component_names.begin(), _component_names.end(), components[c]) != _component_names.end()))
        for (unsigned int j=0; j<_component_names.size(); ++j)
          mask[j] = (_component_names[j] == components[c] || mask[j]);
      else
        {
          try
            {
              unsigned int m = Utilities::string_to_int(components[c]);
              AssertThrow(m < n_components, ExcWrongComponent(m,n_components));
              mask[m] = true;
            }
          catch (std::exception &exc)
            AssertThrow(false, ExcWrongVariable(components[c],_component_names));
        }

    }

  // boundary
  for (unsigned int c=0; c<boundary_components.size(); ++c)
    {

      if ((std::find(_component_names.begin(), _component_names.end(), boundary_components[c]) != _component_names.end()))
        for (unsigned int j=0; j<_component_names.size(); ++j)
          boundary_mask[j] = (_component_names[j] == boundary_components[c] || mask[j]);

      else
        {
          try
            {
              unsigned int m = Utilities::string_to_int(boundary_components[c]);
              AssertThrow(m < n_components, ExcWrongComponent(m,n_components));
              boundary_mask[m] = true;
            }
          catch (std::exception &exc)
            AssertThrow(false, ExcWrongVariable(boundary_components[c],_component_names));
        }

    }

}



template <int dim, int spacedim>
ComponentMask ParsedZeroAverageConstraints<dim,spacedim>::get_mask() const
{
  return ComponentMask(mask);
}


template <int dim, int spacedim>
ComponentMask ParsedZeroAverageConstraints<dim,spacedim>::get_boundary_mask() const
{
  return ComponentMask(boundary_mask);
}


template <int dim, int spacedim>
void ParsedZeroAverageConstraints<dim,spacedim>::declare_parameters(ParameterHandler &prm)
{
  if (str_component_names != "")
    add_parameter(prm, &_component_names, "Known component names", str_component_names,
                  Patterns::List(Patterns::Anything(),1,n_components,","),
                  "These variables can be used to set the corrisponding component mask, "
                  "instead of specifying each component number");

  else
    {
      std::vector<std::string> cn(n_components, "u");
      add_parameter(prm, &_component_names, "Known component names", print(cn),
                    Patterns::List(Patterns::Anything(),1,n_components,","),
                    "These variables can be used to set the corrisponding component mask, "
                    "instead of specifying each component number");

    }
  add_parameter(prm, &components, "Zero average on whole domain", str_components,
                Patterns::List(Patterns::Anything(),0,n_components,","),
                "Pattern to be used: "
                "0,2,p \n"
                "You can specify the components either by numbers "
                "or by the corrisponding variable name, which are parsed at "
                "construction time. ");

  add_parameter(prm, &boundary_components, "Zero average on boundary", str_boundary_components,
                Patterns::List(Patterns::Anything(),0,n_components,","),
                "Pattern to be used: "
                "0,2,p \n"
                "You can specify the components either by numbers "
                "or by the corrisponding variable name, which are parsed at "
                "construction time. ");

}

template <>
void ParsedZeroAverageConstraints<1,1>::
internal_zero_average_constraints(const DoFHandler<1,1> &,
                                  const ComponentMask,
                                  const bool,
                                  ConstraintMatrix &) const
{
  ExcImpossibleInDim(1);
}

template <>
void ParsedZeroAverageConstraints<1,2>::
internal_zero_average_constraints(const DoFHandler<1,2> &,
                                  const ComponentMask,
                                  const bool,
                                  ConstraintMatrix &) const
{
  ExcImpossibleInDim(1);
}

template <>
void ParsedZeroAverageConstraints<1,3>::
internal_zero_average_constraints(const DoFHandler<1,3> &,
                                  const ComponentMask,
                                  const bool,
                                  ConstraintMatrix &) const
{
  ExcImpossibleInDim(1);
}


template <int dim, int spacedim>
void ParsedZeroAverageConstraints<dim,spacedim>::
internal_zero_average_constraints(const DoFHandler<dim,spacedim> &dof_handler,
                                  const ComponentMask mask,
                                  const bool at_boundary,
                                  ConstraintMatrix &constraints) const
{
  std::vector<bool> constrained_dofs (dof_handler.n_dofs(), false);

  if (at_boundary)
    DoFTools::extract_boundary_dofs (dof_handler,
                                     mask,
                                     constrained_dofs);
  else
    DoFTools::extract_dofs (dof_handler,
                            mask,
                            constrained_dofs);


  const unsigned int first_dof
    = std::distance (constrained_dofs.begin(),
                     std::find (constrained_dofs.begin(),
                                constrained_dofs.end(),
                                true) );

  constraints.add_line (first_dof);
  for (unsigned int i=first_dof+1; i<dof_handler.n_dofs(); ++i)
    if (constrained_dofs[i] == true)
      constraints.add_entry(first_dof, i, -1);
}


template <int dim, int spacedim>
void ParsedZeroAverageConstraints<dim,spacedim>::
apply_zero_average_constraints(const DoFHandler<dim,spacedim> &dof_handler,
                               ConstraintMatrix &constraints) const
{
  for (unsigned int i=0; i<n_components; ++i)
    {
      std::vector<bool> m(n_components,false);
      std::vector<bool> m_boundary(n_components,false);

      if (boundary_mask[i])
        {
          m_boundary[i] = true;
          internal_zero_average_constraints(dof_handler,
                                            ComponentMask(m_boundary),
                                            true,
                                            constraints);
        }

      if (mask[i])
        {
          m[i] = true;
          internal_zero_average_constraints(dof_handler,
                                            ComponentMask(m),
                                            false,
                                            constraints);
        }
    }
}


D2K_NAMESPACE_CLOSE

template class deal2lkit::ParsedZeroAverageConstraints<1,1>;
template class deal2lkit::ParsedZeroAverageConstraints<1,2>;
template class deal2lkit::ParsedZeroAverageConstraints<1,3>;
template class deal2lkit::ParsedZeroAverageConstraints<2,2>;
template class deal2lkit::ParsedZeroAverageConstraints<2,3>;
template class deal2lkit::ParsedZeroAverageConstraints<3,3>;
