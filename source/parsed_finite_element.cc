#include <deal2lkit/parsed_finite_element.h>
#include <deal2lkit/utilities.h>
#include <deal.II/fe/fe_tools.h>

template <int dim, int spacedim>
ParsedFiniteElement<dim, spacedim>::ParsedFiniteElement(const std::string &name,
                                                        const std::string &default_name,
                                                        const std::string &default_component_names,
                                                        const unsigned int n_components,
                                                        const std::string &default_coupling,
                                                        const std::string &default_preconditioner_coupling) :
  ParameterAcceptor(name),
  _n_components(n_components),
  fe_name(default_name),
  default_component_names(default_component_names),
  default_coupling(default_coupling),
  default_preconditioner_coupling(default_preconditioner_coupling)
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

  add_parameter(prm, &coupling_int,
                "Block coupling", default_coupling,
                Patterns::List(Patterns::List(Patterns::Integer(0,3),0,
                                              (_n_components ? _n_components: numbers::invalid_unsigned_int), ","),
                               0, (_n_components ? _n_components: numbers::invalid_unsigned_int), ";"),
                "Coupling between the blocks of the finite elements in the system:\n"
                " 0: No coupling\n"
                " 1: Full coupling\n"
                " 2: Coupling only on faces\n");


  add_parameter(prm, &preconditioner_coupling_int,
                "Preconditioner block coupling", default_preconditioner_coupling,
                Patterns::List(Patterns::List(Patterns::Integer(0,3),0,
                                              (_n_components ? _n_components: numbers::invalid_unsigned_int), ","),
                               0, (_n_components ? _n_components: numbers::invalid_unsigned_int), ";"),
                "Coupling between the blocks of the finite elements in the preconditioner:\n"
                " 0: No coupling\n"
                " 1: Full coupling\n"
                " 2: Coupling only on faces\n");
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

  coupling = to_coupling(coupling_int);
  preconditioner_coupling = to_coupling(preconditioner_coupling_int);
}



template<int dim, int spacedim>
Table<2,DoFTools::Coupling> ParsedFiniteElement<dim,spacedim>::to_coupling(const std::vector<std::vector<unsigned int> > &coupling_table) const
{
  const unsigned int nc = n_components();
  const unsigned int nb = n_blocks();

  Table<2,DoFTools::Coupling> out_coupling(nc, nc);

  std::vector<DoFTools::Coupling> m(3);
  m[0] = DoFTools::none;
  m[1] = DoFTools::always;
  m[2] = DoFTools::nonzero;

  if (coupling_table.size() == nc)
    for (unsigned int i=0; i<nc; ++i)
      {
        AssertThrow(coupling_table[i].size() == nc, ExcDimensionMismatch(coupling_table[i].size(), nc));
        for (unsigned int j=0; j<nc; ++j)
          out_coupling[i][j] = m[coupling_table[i][j]];
      }
  else if (coupling_table.size() == nb)
    for (unsigned int i=0; i<nc; ++i)
      {
        AssertThrow(coupling_table[component_blocks[i]].size() == nb,
                    ExcDimensionMismatch(coupling_table[component_blocks[i]].size(), nb));
        for (unsigned int j=0; j<nc; ++j)
          out_coupling[i][j] = m[coupling_table[component_blocks[i]][component_blocks[j]]];
      }
  else if (coupling_table.size() == 0)
    for (unsigned int i=0; i<nc; ++i)
      {
        for (unsigned int j=0; j<nc; ++j)
          out_coupling[i][j] = m[1];
      }
  else
    AssertThrow(false, ExcMessage("You tried to construct a coupling with the wrong number of elements."));

  return out_coupling;
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
const Table<2, DoFTools::Coupling> &ParsedFiniteElement<dim,spacedim>::get_coupling() const
{
  return coupling;
}

template<int dim, int spacedim>
const Table<2, DoFTools::Coupling> &ParsedFiniteElement<dim,spacedim>::get_preconditioner_coupling() const
{
  return preconditioner_coupling;
}


template class ParsedFiniteElement<1,1>;
template class ParsedFiniteElement<1,2>;
template class ParsedFiniteElement<1,3>;
template class ParsedFiniteElement<2,2>;
template class ParsedFiniteElement<2,3>;
template class ParsedFiniteElement<3,3>;

