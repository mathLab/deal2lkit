#ifndef __dealii_sak_parsed_mapped_function_h
#define __dealii_sak_parsed_mapped_function_h

#include <deal.II/base/exceptions.h>
#include "parameter_acceptor.h"
#include "parsed_function.h"
#include <deal.II/fe/component_mask.h>
#include <algorithm>
#include <map>

using namespace dealii;

/**
 * ParsedMappedFunctions object.
 * It allows you to set a mapped functions, i.e., parsed functions
 * acting on specified id (boundary_id, material_id,etc..) and on
 * specified components.
 *
 * Dirichlet Boundary conditions or Neumann boundary conditions
 *
 * can be set using this class.
 */
template <int spacedim, int n_components>
class ParsedMappedFunctions : public ParameterAcceptor
{
public:
//  ParsedMappedFunctions  () {};
  ParsedMappedFunctions  (const std::string &name = "",
                          const std::string &component_names = "u",
                          const std::string &default_id_components = "", // 0: 0; 1; 2 , 1: 3
                          const std::string &default_id_functions = "",  // 0: x; y; z, 1: sin(k)
                          const std::string &default_constants = "");

  shared_ptr<ParsedFunction<spacedim,n_components> > get_mapped_function (const unsigned int &id) const;

  ComponentMask get_mapped_mask (const unsigned int &id) const;

  std::vector<unsigned int> get_mapped_ids() const;


  virtual void declare_parameters (ParameterHandler &prm);

  /// An entry with this id does not exist in this object.
  DeclException1(ExcIdNotFound, unsigned int,
                 << "No entry with the id " << arg1 << " exists.");
private:

  void split_id_components(const std::string &parsed_id_components);
  void split_id_functions(const std::string &parsed_id_functions,
                          const std::string &constants);


  std::string name;
  std::vector<std::string> component_names;
  std::vector<unsigned int> id_defined_components;
  std::vector<unsigned int> id_defined_functions;
  std::map<unsigned int, ComponentMask> id_components;
  std::map<unsigned int, shared_ptr<ParsedFunction<spacedim,n_components> > > id_functions;
  std::map<unsigned int, std::pair<ComponentMask, shared_ptr<ParsedFunction<spacedim,n_components> > > > mapped_functions;

};

template <int spacedim, int n_components>
ParsedMappedFunctions<spacedim,n_components>::ParsedMappedFunctions(const std::string &parsed_name,
    const std::string &parsed_component_names,
    const std::string &parsed_id_components,
    const std::string &parsed_id_functions,
    const std::string &parsed_constants):
  ParameterAcceptor(parsed_name),
  name (parsed_name)
{
  component_names = Utilities::split_string_list(parsed_component_names, ',');
  split_id_components(parsed_id_components);
  split_id_functions(parsed_id_functions,parsed_constants);

  AssertDimension(id_components.size(), id_functions.size());

  AssertDimension(id_defined_components.size(), id_defined_functions.size());

  Assert((std::is_permutation(id_defined_components.begin(),
                              id_defined_components.end(),
                              id_defined_functions.begin())),
         ExcMessage("Ids associated to components and to functions are not the same."));

  typedef std::map<unsigned int, ComponentMask>::iterator it_type;
  for (it_type it=id_components.begin(); it != id_components.end(); ++it)
    {
      std::pair<ComponentMask, shared_ptr<ParsedFunction<spacedim,n_components> > > mapped(it->second, id_functions[it->first]);

      mapped_functions[it->first] = mapped;

    }

};

template <int spacedim, int n_components>
void ParsedMappedFunctions<spacedim,n_components>::split_id_components(const std::string &parsed_id_components)
{
  std::vector<std::string> idcomponents;

  idcomponents = Utilities::split_string_list(parsed_id_components, ',');

  for (unsigned int i=0; i<idcomponents.size(); ++i)
    {
      std::vector<std::string> id_comp;
      std::vector<std::string> components;
      std::vector<bool> mask(n_components,false);

      id_comp = Utilities::split_string_list(idcomponents[i], ':');

      unsigned int id = Utilities::string_to_int(id_comp[0]);
      id_defined_components.push_back(id);

      components = Utilities::split_string_list(id_comp[1], ';');
      for (unsigned int c=0; c<components.size(); ++c)
        {
          unsigned int m = Utilities::string_to_int(components[c]);
          mask[m] = true;
        }
      id_components[id] = ComponentMask(mask);
    }
}


template <int spacedim, int n_components>
void ParsedMappedFunctions<spacedim,n_components>::split_id_functions(const std::string &parsed_id_functions,
    const std::string &constants)
{
  std::vector<std::string> idfunctions;

  idfunctions = Utilities::split_string_list(parsed_id_functions, ',');

  for (unsigned int i=0; i<idfunctions.size(); ++i)
    {
      std::vector<std::string> id_func;

      id_func = Utilities::split_string_list(idfunctions[i], ':');

      unsigned int id = Utilities::string_to_int(id_func[0]);
      id_defined_functions.push_back(id);

      std::string function_name = name + " acting on id " + Utilities::int_to_string(id);
      shared_ptr<ParsedFunction<spacedim,n_components> > pf;
      pf = SP(new ParsedFunction<spacedim,n_components>(function_name, id_func[1], constants));

      id_functions[id] = pf;
    }
}

template <int spacedim, int n_components>
shared_ptr<ParsedFunction<spacedim,n_components> > ParsedMappedFunctions<spacedim,n_components>::get_mapped_function(const unsigned int &id) const
{
  Assert( mapped_functions.find(id) != mapped_functions.end(),
          ExcIdNotFound(id));
  return mapped_functions.at(id).second;
}

template <int spacedim, int n_components>
ComponentMask ParsedMappedFunctions<spacedim,n_components>::get_mapped_mask(const unsigned int &id) const
{
  Assert( mapped_functions.find(id) != mapped_functions.end(),
          ExcIdNotFound(id));
  return mapped_functions.at(id).first;
}

template <int spacedim, int n_components>
std::vector<unsigned int> ParsedMappedFunctions<spacedim,n_components>::get_mapped_ids() const
{
  return id_defined_components;
}

template <int spacedim, int n_components>
void ParsedMappedFunctions<spacedim,n_components>::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &name,
                "Mapped Functions", name, Patterns::Anything(),
                "Name of the mapped functions");

}

#endif
