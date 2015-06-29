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
 * Dirichlet Boundary conditions or Neumann boundary conditions
 *
 * can be set using this class.
 */
template <int spacedim, int n_components>
class ParsedMappedFunctions : public ParameterAcceptor
{
public:
//  ParsedMappedFunctions  () {};
  ParsedMappedFunctions  (const std::string &name = "Mapped Functions",
                          const std::string &component_names = "",
                          const std::string &default_id_components = "0=ALL",
                          const std::string &default_id_functions = "",
                          const std::string &default_constants = "");

  shared_ptr<dealii::Functions::ParsedFunction<spacedim> > get_mapped_function (const unsigned int &id) const;

  ComponentMask get_mapped_mask (const unsigned int &id) const;

  std::vector<unsigned int> get_mapped_ids() const;


  virtual void declare_parameters (ParameterHandler &prm);
  virtual void parse_parameters_call_back ();

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

private:

  void split_id_components(const std::string &parsed_idcomponents);
  void split_id_functions(const std::string &parsed_idfunctions,
                          const std::string &constants);


  std::string name;
  std::string str_id_components;
  std::string str_id_functions;
  std::string str_component_names;
  std::string str_constants;
  std::vector<std::string> _component_names;
  std::vector<unsigned int> ids;
  std::map<unsigned int, ComponentMask> id_components;
  std::map<unsigned int, shared_ptr<dealii::Functions::ParsedFunction<spacedim> > > id_functions;
  std::map<unsigned int, std::pair<ComponentMask, shared_ptr<dealii::Functions::ParsedFunction<spacedim> > > > mapped_functions;

};

template <int spacedim, int n_components>
ParsedMappedFunctions<spacedim,n_components>::ParsedMappedFunctions(const std::string &parsed_name,
    const std::string &parsed_component_names,
    const std::string &parsed_id_components,
    const std::string &parsed_id_functions,
    const std::string &parsed_constants):
  ParameterAcceptor(parsed_name),
  name (parsed_name),
  str_component_names (parsed_component_names),
  str_id_components (parsed_id_components),
  str_id_functions (parsed_id_functions),
  str_constants (parsed_constants)
{};

template <int spacedim, int n_components>
void ParsedMappedFunctions<spacedim,n_components>::parse_parameters_call_back()
{
  split_id_components(str_id_components);
  split_id_functions(str_id_functions,str_constants);
  AssertDimension(id_components.size(), id_functions.size());

  typedef std::map<unsigned int, ComponentMask>::iterator it_type;
  for (it_type it=id_components.begin(); it != id_components.end(); ++it)
    {
      std::pair<ComponentMask, shared_ptr<dealii::Functions::ParsedFunction<spacedim> > > mapped(it->second, id_functions[it->first]);

      mapped_functions[it->first] = mapped;

    }
}

template <int spacedim, int n_components>
void ParsedMappedFunctions<spacedim,n_components>::split_id_components(const std::string &parsed_idcomponents)
{
  std::vector<std::string> idcomponents;

  idcomponents = Utilities::split_string_list(parsed_idcomponents, '%');

  for (unsigned int i=0; i<idcomponents.size(); ++i)
    {
      std::vector<std::string> id_comp;
      std::vector<std::string> components;
      std::vector<bool> mask(n_components,false);

      id_comp = Utilities::split_string_list(idcomponents[i], '=');

      unsigned int id = Utilities::string_to_int(id_comp[0]);
      ids.push_back(id);

      components = Utilities::split_string_list(id_comp[1], ';');
      Assert(components.size() <= n_components,
             ExcMessage("Wrong number of components specified for id " + id_comp[0]));
      for (unsigned int c=0; c<components.size(); ++c)
        {
          if (components[c] == "ALL")
            {
              for (unsigned int j=0; j<n_components; ++j)
                mask[j] = true;
              break;
            }

          Assert(_component_names.size() == n_components,
                 ExcMessage("Parsed components name must match the number of components of the problem"));

          if ((std::find(_component_names.begin(), _component_names.end(), components[c]) != _component_names.end()))
            for (unsigned int j=0; j<_component_names.size(); ++j)
              mask[j] = (_component_names[j] == components[c] || mask[j]);
          else
            {
              try
                {
                  unsigned int m = Utilities::string_to_int(components[c]);
                  Assert(m < n_components, ExcWrongComponent(id,m,n_components));
                  mask[m] = true;
                }
              catch (std::exception &exc)
                Assert(false, ExcWrongVariable(id,components[c],_component_names));
            }
        }
      id_components[id] = ComponentMask(mask);
    }
}

template <int spacedim, int n_components>
void ParsedMappedFunctions<spacedim,n_components>::split_id_functions(const std::string &parsed_idfunctions,
    const std::string &constants)
{
  std::vector<unsigned int> id_defined_functions;

  // if it is empty a ZeroFunction<dim>(n_components) is applied on the
  // parsed ids in the components
  if (parsed_idfunctions == "")
    {
      for (unsigned int i=0; i<ids.size(); ++i)
        {
          id_defined_functions.push_back(ids[i]);
          shared_ptr<dealii::Functions::ParsedFunction<spacedim> > ptr;

          ParameterHandler internal_prm;
          dealii::Functions::ParsedFunction<spacedim>::declare_parameters(internal_prm, n_components);
          ptr = SP(new dealii::Functions::ParsedFunction<spacedim>(n_components));
          ptr->parse_parameters(internal_prm);

          id_functions[ids[i]] = ptr;
        }
    }
  else
    {


      std::vector<std::string> idfunctions;
      idfunctions = Utilities::split_string_list(parsed_idfunctions, '%');

      for (unsigned int i=0; i<idfunctions.size(); ++i)
        {
          std::vector<std::string> id_func;

          id_func = Utilities::split_string_list(idfunctions[i], '=');

          unsigned int id = Utilities::string_to_int(id_func[0]);

          // check if the current id is also defined in id_components
          Assert((std::find(ids.begin(), ids.end(), id) != ids.end()),
                 ExcIdNotMatch(id));
          id_defined_functions.push_back(id);

          shared_ptr<dealii::Functions::ParsedFunction<spacedim> > ptr;

          ParameterHandler internal_prm;
          dealii::Functions::ParsedFunction<spacedim>::declare_parameters(internal_prm, n_components);
          internal_prm.set("Function expression", id_func[1]);
          internal_prm.set("Function constants", constants);
          ptr = SP(new dealii::Functions::ParsedFunction<spacedim>(n_components));
          ptr->parse_parameters(internal_prm);

          id_functions[id] = ptr;
        }
    }

  // check if the number of ids defined in id_components and id_functions are the same
  Assert(ids.size() == id_defined_functions.size(),
         ExcMessage("Ids associated to components and to functions are not the same."));

}

template <int spacedim, int n_components>
shared_ptr<dealii::Functions::ParsedFunction<spacedim> > ParsedMappedFunctions<spacedim,n_components>::get_mapped_function(const unsigned int &id) const
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
  return ids;
}

template <int spacedim, int n_components>
void ParsedMappedFunctions<spacedim,n_components>::declare_parameters(ParameterHandler &prm)
{
  if (str_component_names != "")
    add_parameter(prm, &_component_names, "Known component names", str_component_names,
                  Patterns::List(Patterns::Anything(),1,n_components,","),
                  "These variables can be used to set the corrisponding component mask,"
                  "instead of specifying each component number");

  else
    {
      std::vector<std::string> cn(n_components, "u");
      add_parameter(prm, &_component_names, "Known component names", print(cn),
                    Patterns::List(Patterns::Anything(),1,n_components,","),
                    "These variables can be used to set the corrisponding component mask,"
                    "instead of specifying each component number");

    }
  add_parameter(prm, &str_id_components, "IDs and component masks", str_id_components,
                Patterns::Anything(),
                "Pattern to be used"
                "id followed by '=' component masks separated by ';'"
                "each couple of id and mask is separated by '%'"
                "0=0;1;2 % 4=u;p % 2=3 % 5=ALL"
                "You can specify the components either by numbers "
                "or by the corrisponding variable name, which are parsed at"
                "construction time."
                "The keyword 'ALL' means all the components.");

  add_parameter(prm, &str_id_functions, "IDs and expressions", str_id_functions,
                Patterns::Anything(),
                "Pattern to be used"
                "id followed by '=' component separated by ';'"
                "each couple of id and expression _id_functions separated by '%'"
                "0=x;y;k;0 % 4=sin(x);cos(y);2*k;1 % 2=0;0;0;0"
                "If it is left empty, a ZeroFunction<dim>(n_components)"
                "is applied on the parsed ids in the components. ");

  add_parameter(prm, &str_constants , "Used constants", str_constants, Patterns::Anything(),
                "Costants which are employed in the definitions of the function expressions."
                "The pattern to be used is"
                "constant_name=value , other_constant=other_value");

}

#endif
