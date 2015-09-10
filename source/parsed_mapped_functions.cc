#include "parsed_mapped_functions.h"

template <int spacedim, int n_components>
ParsedMappedFunctions<spacedim,n_components>::ParsedMappedFunctions(const std::string &parsed_name,
    const std::string &parsed_component_names,
    const std::string &parsed_id_components,
    const std::string &parsed_id_functions,
    const std::string &parsed_constants):
  ParameterAcceptor(parsed_name),
  name (parsed_name),
  str_id_components (parsed_id_components),
  str_id_functions (parsed_id_functions),
  str_component_names (parsed_component_names),
  str_constants (parsed_constants)
{}

template <int spacedim, int n_components>
void ParsedMappedFunctions<spacedim,n_components>::add_normal_components()
{
  std::vector<std::string> var = unique(_component_names);
  for (unsigned int i=0; i<var.size(); ++i)
    {
      int n = std::count(_component_names.begin(),_component_names.end(), var[i]);
      if (n>1)
        {
          _normal_components.push_back(var[i]+".N");
          int pos = std::find(_component_names.begin(),_component_names.end(), var[i]) - _component_names.begin();
          std::pair<unsigned int, std::string> nc (pos, var[i]+".N");
          normal_components.push_back(nc);
        }
    }
}

template <int spacedim, int n_components>
void ParsedMappedFunctions<spacedim,n_components>::parse_parameters_call_back()
{
  add_normal_components();
  _all_components = _component_names;
  for (unsigned int i=0; i<_normal_components.size(); ++i)
    {
      _all_components.push_back(_normal_components[i]);
    }
  split_id_components(str_id_components);
  split_id_functions(str_id_functions,str_constants);
  AssertThrow(id_components.size() == id_functions.size(),
              ExcIdsMismatch(id_components.size(), id_functions.size()));

  typedef std::map<unsigned int, ComponentMask>::iterator it_type;
  for (it_type it=id_components.begin(); it != id_components.end(); ++it)
    {
      std::pair<ComponentMask, shared_ptr<dealii::Functions::ParsedFunction<spacedim> > > mapped(it->second, id_functions[it->first]);

      mapped_functions[it->first] = mapped;

    }

  set_normal_functions();

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
      AssertThrow(components.size() <= n_components,
                  ExcMessage("Wrong number of components specified for id " + id_comp[0]));
      for (unsigned int c=0; c<components.size(); ++c)
        {
          if (components[c] == "ALL")
            {
              for (unsigned int j=0; j<n_components; ++j)
                mask[j] = true;
              break;
            }

          AssertThrow(_component_names.size() == n_components,
                      ExcMessage("Parsed components name must match the number of components of the problem"));

          if ((std::find(_component_names.begin(), _component_names.end(), components[c]) != _component_names.end()))
            for (unsigned int j=0; j<_component_names.size(); ++j)
              mask[j] = (_component_names[j] == components[c] || mask[j]);

          else if ((std::find(_normal_components.begin(), _normal_components.end(), components[c]) != _normal_components.end()))
            {
              normal_ids.push_back(id);
              for (unsigned int j=0; j<normal_components.size(); ++j)
                if (normal_components[j].second == components[c])
                  {
                    mapped_normal_components[components[c]].second = normal_components[j].first;
                    mapped_normal_components[components[c]].first.push_back(id);
                  }
            }
          else
            {
              try
                {
                  unsigned int m = Utilities::string_to_int(components[c]);
                  AssertThrow(m < n_components, ExcWrongComponent(id,m,n_components));
                  mask[m] = true;
                }
              catch (std::exception &exc)
                AssertThrow(false, ExcWrongVariable(id,components[c],_all_components));
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
          std::string str;
          for (unsigned int i=0; i<n_components-1; ++i)
            str += "0;";
          str+="0";

          id_str_functions[ids[i]] = str;

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
          id_str_functions[id] = id_func[1];

          // check if the current id is also defined in id_components
          AssertThrow((std::find(ids.begin(), ids.end(), id) != ids.end()),
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
  AssertThrow(ids.size() == id_defined_functions.size(),
              ExcMessage("Ids associated to components and to functions are not the same."));

}

template <int spacedim, int n_components>
shared_ptr<dealii::Functions::ParsedFunction<spacedim> >
ParsedMappedFunctions<spacedim,n_components>::get_mapped_normal_function(const unsigned int &id, const unsigned int &fcv) const
{
  std::pair<unsigned int, unsigned int> id_fcv (id,fcv);
  Assert( _normal_functions.find(id_fcv) != _normal_functions.end(),
          ExcIdNotFound(id));

  return _normal_functions.at(id_fcv);
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
std::vector<unsigned int> ParsedMappedFunctions<spacedim,n_components>::get_mapped_normal_ids() const
{
  return normal_ids;
}

template <int spacedim, int n_components>
void ParsedMappedFunctions<spacedim,n_components>::declare_parameters(ParameterHandler &prm)
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
  add_parameter(prm, &str_id_components, "IDs and component masks", str_id_components,
                Patterns::Anything(),
                "Pattern to be used "
                "id followed by '=' component masks separated by ';' "
                "each couple of id and mask is separated by '%' "
                "0=0;1;2 % 4=u;p % 2=3 % 5=ALL "
                "You can specify the components either by numbers "
                "or by the corrisponding variable name, which are parsed at "
                "construction time. "
                "The keyword 'ALL' means all the components. "
                "Normal component is referred with suffix N "
                "e.g. uN means the normal component of u. "
                "note that the normal component can be set only "
                "for a vector variable.");

  add_parameter(prm, &str_id_functions, "IDs and expressions", str_id_functions,
                Patterns::Anything(),
                "Pattern to be used  "
                "id followed by '=' component separated by ';' "
                "each couple of id and expression _id_functions separated by '%' "
                "0=x;y;k;0 % 4=sin(x);cos(y);2*k;1 % 2=0;0;0;0 "
                "If it is left empty, a ZeroFunction<dim>(n_components) "
                "is applied on the parsed ids in the components. ");

  add_parameter(prm, &str_constants , "Used constants", str_constants, Patterns::Anything(),
                "Costants which are employed in the definitions of the function expressions. "
                "The pattern to be used is "
                "constant_name=value , other_constant=other_value");

}

template <int spacedim, int n_components>
bool ParsedMappedFunctions<spacedim,n_components>::acts_on_id(unsigned int &id) const
{
  return id_components.find(id) != id_components.end();
}

template <int spacedim, int n_components>
void ParsedMappedFunctions<spacedim,n_components>::set_time(const double &t)
{
  typedef typename std::map<unsigned int, shared_ptr<dealii::Functions::ParsedFunction<spacedim> > >::iterator it_type;
  for (it_type it=id_functions.begin(); it != id_functions.end(); ++it)
    it->second->set_time(t);
}

template <int spacedim, int n_components>
void ParsedMappedFunctions<spacedim,n_components>::set_normal_functions()
{
  typedef std::map<std::string, std::pair<std::vector<unsigned int>, unsigned int > >::iterator it_type;

  for (it_type it=mapped_normal_components.begin(); it != mapped_normal_components.end(); ++it)
    {

      std::vector<unsigned int> normal_ids = (it->second).first;
      unsigned int fcv = (it->second).second; // unsigned int first component vector

      for (unsigned int i=0; i<normal_ids.size(); ++i)
        {
          std::vector<std::string> normal_func;
          normal_func = Utilities::split_string_list(id_str_functions[normal_ids[i]], ';');
          shared_ptr<dealii::Functions::ParsedFunction<spacedim> > normal_ptr;

          ParameterHandler normal_prm;
          dealii::Functions::ParsedFunction<spacedim>::declare_parameters(normal_prm, spacedim);
          std::string str_normal_func;
          Assert(spacedim>1, ExcNotImplemented());
          for (unsigned int i=0; i<spacedim-1; ++i)
            str_normal_func += normal_func[fcv+i] + ";";
          str_normal_func += normal_func[fcv+spacedim-1];

          normal_prm.set("Function expression", str_normal_func);
          normal_prm.set("Function constants", str_constants);
          normal_ptr = SP(new dealii::Functions::ParsedFunction<spacedim>(spacedim));
          normal_ptr->parse_parameters(normal_prm);
          std::pair<unsigned int, unsigned int> id_fcv (normal_ids[i],fcv);
          _normal_functions[id_fcv]=normal_ptr;
        }
    }
}

template class ParsedMappedFunctions<1,1>;

template class ParsedMappedFunctions<2,1>;
template class ParsedMappedFunctions<2,2>;
template class ParsedMappedFunctions<2,3>;
template class ParsedMappedFunctions<2,4>;
template class ParsedMappedFunctions<2,5>;
template class ParsedMappedFunctions<2,6>;
template class ParsedMappedFunctions<2,7>;
template class ParsedMappedFunctions<2,8>;

template class ParsedMappedFunctions<3,1>;
template class ParsedMappedFunctions<3,2>;
template class ParsedMappedFunctions<3,3>;
template class ParsedMappedFunctions<3,4>;
template class ParsedMappedFunctions<3,5>;
template class ParsedMappedFunctions<3,6>;
template class ParsedMappedFunctions<3,7>;
template class ParsedMappedFunctions<3,8>;
