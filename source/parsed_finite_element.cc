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

#include <deal2lkit/parsed_finite_element.h>
#include <deal2lkit/utilities.h>
#include <deal.II/fe/fe_tools.h>

#include <algorithm> // std::find

D2K_NAMESPACE_OPEN

template <int dim, int spacedim>
ParsedFiniteElement<dim, spacedim>::ParsedFiniteElement(const std::string &name,
                                                        const std::string &default_name,
                                                        const std::string &default_component_names,
                                                        const unsigned int n_components) :
  ParameterAcceptor(name),
  _n_components(n_components),
  fe_name(default_name),
  default_component_names(default_component_names)
{
  component_names = Utilities::split_string_list(default_component_names);
  parse_parameters_call_back();
}

template <int dim, int spacedim>
void ParsedFiniteElement<dim, spacedim>::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &fe_name,
                "Finite element space", fe_name,
                Patterns::Anything(),
                "The finite element space to use. For vector "
                "finite elements use the notation "
                "FESystem[FE_Q(2)^2-FE_DGP(1)] (e.g. Navier-Stokes). ");

  add_parameter(prm, &component_names,
                "Blocking of the finite element", default_component_names,
                // This ensures that an assert is thrown if you try to
                // read something with the wrong number of components
                Patterns::List(Patterns::Anything(),
                               (_n_components ? _n_components: 1),
                               (_n_components ? _n_components: numbers::invalid_unsigned_int)),
                "How to partition the finite element. This information can be used "
                "to construct block matrices and vectors, as well as to create "
                "names for solution vectors, or error tables. A repeated component "
                "is interpreted as a vector field, with dimension equal to the "
                "number of repetitions (up to 3). This is used in conjunction "
                "with a ParsedFiniteElement class, to generate arbitrary "
                "finite dimensional spaces.");
}

template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
ParsedFiniteElement<dim, spacedim>::operator()() const
{
  return FETools::get_fe_by_name<dim,spacedim>(fe_name);
}


template<int dim, int spacedim>
void ParsedFiniteElement<dim,spacedim>::parse_parameters_call_back()
{
  component_blocks.resize(component_names.size());
  block_names.resize(component_names.size());
  unsigned int j=0;
  for (unsigned int i=0; i<component_names.size(); ++i)
    {
      if ( (i>0) && (component_names[i-1] != component_names[i]) )
        j++;
      component_blocks[i] = j;
      block_names[j] = component_names[i];
    }
  block_names.resize(j+1);
  FiniteElement<dim,spacedim> *fe = (*this)();
  const unsigned int nc = fe->n_components();
  delete fe;
  AssertThrow(component_names.size() == nc,
              ExcInternalError("Generated FE has the wrong number of components."));
}


template<int dim, int spacedim>
unsigned int ParsedFiniteElement<dim,spacedim>::n_components() const
{
  return component_names.size();
}


template<int dim, int spacedim>
unsigned int ParsedFiniteElement<dim,spacedim>::n_blocks() const
{
  return block_names.size();
}


template<int dim, int spacedim>
std::string ParsedFiniteElement<dim,spacedim>::get_component_names() const
{
  return print(component_names);
}


template<int dim, int spacedim>
std::string ParsedFiniteElement<dim,spacedim>::get_block_names() const
{
  return print(block_names);
}


template<int dim, int spacedim>
std::vector<unsigned int> ParsedFiniteElement<dim,spacedim>::get_component_blocks() const
{
  return component_blocks;
}

template<int dim, int spacedim>
unsigned int ParsedFiniteElement<dim,spacedim>::get_component_position(const std::string &var) const
{
  unsigned int pos_counter = 0;
  auto pos_it = std::find (component_names.begin(), component_names.end(), var);
  Assert(pos_it != component_names.end(),
         ExcInternalError("Component not found!"));
  return pos_it - component_names.begin();
}

template<int dim, int spacedim>
bool ParsedFiniteElement<dim,spacedim>::is_vectorial(const std::string &var) const
{
  auto pos_it = std::find (component_names.begin(), component_names.end(), var);
  Assert(pos_it != component_names.end(),
         ExcInternalError("Component not found!"));
  pos_it++;
  return (*pos_it == var);
}

D2K_NAMESPACE_CLOSE


template class deal2lkit::ParsedFiniteElement<1,1>;
template class deal2lkit::ParsedFiniteElement<1,2>;
template class deal2lkit::ParsedFiniteElement<1,3>;
template class deal2lkit::ParsedFiniteElement<2,2>;
template class deal2lkit::ParsedFiniteElement<2,3>;
template class deal2lkit::ParsedFiniteElement<3,3>;

