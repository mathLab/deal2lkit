#include "parsed_finite_element.h"
#include <deal.II/fe/fe_tools.h>

namespace
{
  std::string join(const std::vector<std::string> &list)
  {
    std::string ret = (list.size() ? list[0]: "");
    for (unsigned int i=1; i<list.size(); ++i)
      ret += ", " + list[i];
    return ret;
  }
}


template <int dim, int spacedim>
ParsedFiniteElement<dim, spacedim>::ParsedFiniteElement(const std::string &name,
                                                        const std::string &default_name,
                                                        const std::string &default_component_names,
                                                        const unsigned int n_components) :
  ParameterAcceptor(name),
  _n_components(n_components),
  fe_name(default_name),
  default_component_names(default_component_names)
{}

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
ParsedFiniteElement<dim, spacedim>::operator()()
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
  unsigned int nc = fe->n_components();
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
  return join(component_names);
}


template<int dim, int spacedim>
std::string ParsedFiniteElement<dim,spacedim>::get_block_names() const
{
  return join(block_names);
}


template<int dim, int spacedim>
std::vector<unsigned int> ParsedFiniteElement<dim,spacedim>::get_component_blocks() const
{
  return component_blocks;
}


template class ParsedFiniteElement<1,1>;
template class ParsedFiniteElement<1,2>;
template class ParsedFiniteElement<1,3>;
template class ParsedFiniteElement<2,2>;
template class ParsedFiniteElement<2,3>;
template class ParsedFiniteElement<3,3>;

