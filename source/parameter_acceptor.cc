//-----------------------------------------------------------
//
//    Copyright (C) 2015 - 2016 by the deal2lkit authors
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

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/utilities.h>
#include <deal.II/base/point.h>
#include <deal.II/base/revision.h>
#include <deal.II/base/path_search.h>
#include <deal2lkit/revision.h>
#include <fstream>

D2K_NAMESPACE_OPEN


// Static empty class list
std::vector<SmartPointer<ParameterAcceptor> > ParameterAcceptor::class_list;
// Static parameter handler
ParameterHandler ParameterAcceptor::prm;

ParameterAcceptor::ParameterAcceptor(const std::string name) :
  acceptor_id(class_list.size()),
  section_name(name)
{
  SmartPointer<ParameterAcceptor> pt(this, type(*this).c_str());
  class_list.push_back(pt);
}


ParameterAcceptor::~ParameterAcceptor()
{
  class_list[acceptor_id] = 0;
}

std::string ParameterAcceptor::get_section_name() const
{
  return (section_name != "" ? section_name : type(*this));
}


void
ParameterAcceptor::initialize(const std::string filename,
                              const std::string out_filename)
{
  // prm.clear();
  declare_all_parameters(prm);
  if (filename != "")
    {
      // check the extension of input file
      if (filename.substr(filename.find_last_of(".") + 1) == "prm")
        {
          try
            {
              prm.parse_input(filename);
            }
          catch (dealii::PathSearch::ExcFileNotFound)
            {
              std::ofstream out(filename);
              Assert(out, ExcIO());
              prm.print_parameters(out, ParameterHandler::Text);
              out.close();
              prm.parse_input(filename);
            }
        }
      else if (filename.substr(filename.find_last_of(".") + 1) == "xml")
        {
          std::ifstream is(filename);
          if (!is)
            {
              std::ofstream out(filename);
              Assert(out, ExcIO());
              prm.print_parameters(out, ParameterHandler::Text);
              out.close();
              is.clear();
              is.open(filename);
              Assert(is, ExcIO());
            }
          prm.parse_input_from_xml(is);
        }
      else
        AssertThrow(false, ExcMessage("Invalid extension of parameter file. Please use .prm or .xml"));
    }
  parse_all_parameters(prm);
  if (out_filename != "")
    {
      std::ofstream outfile(out_filename.c_str());
      Assert(outfile, ExcIO());
      std::string extension = out_filename.substr(out_filename.find_last_of(".") + 1);

      if ( extension == "prm")
        {
          outfile << "# Parameter file generated with " << std::endl
                  << "# D2K_GIT_BRANCH=       " << D2K_GIT_BRANCH << std::endl
                  << "# D2K_GIT_SHORTREV=     " << D2K_GIT_SHORTREV << std::endl
                  << "# DEAL_II_GIT_BRANCH=   " << DEAL_II_GIT_BRANCH  << std::endl
                  << "# DEAL_II_GIT_SHORTREV= " << DEAL_II_GIT_SHORTREV << std::endl;
          prm.print_parameters(outfile, ParameterHandler::ShortText);
        }
      else if (extension == "xml")
        prm.print_parameters(outfile, ParameterHandler::XML);
      else if (extension == "latex" || extension == "tex")
        prm.print_parameters(outfile, ParameterHandler::LaTeX);
      else
        AssertThrow(false,ExcNotImplemented());
    }
}

void
ParameterAcceptor::clear()
{
  for (unsigned int i=0; i<class_list.size(); ++i)
    if (class_list[i] != NULL)
      {
        class_list[i]->parameters.clear();
        class_list[i]->patterns.clear();
      }
  class_list.clear();
  prm.clear();
}

void
ParameterAcceptor::log_info()
{
  deallog.push("ParameterAcceptor");
  for (unsigned int i=0; i<class_list.size(); ++i)
    {
      deallog << "Class " << i << ":";
      if (class_list[i])
        deallog << class_list[i]->get_section_name() << std::endl;
      else
        deallog << " NULL" << std::endl;
    }
  deallog.pop();
}

void ParameterAcceptor::parse_all_parameters(ParameterHandler &prm)
{
  for (unsigned int i=0; i< class_list.size(); ++i)
    if (class_list[i] != NULL)
      {
        class_list[i]->enter_my_subsection(prm);
        class_list[i]->parse_parameters(prm);
        class_list[i]->parse_parameters_call_back();
        class_list[i]->leave_my_subsection(prm);
      }
}

void ParameterAcceptor::declare_all_parameters(ParameterHandler &prm)
{
  for (unsigned int i=0; i< class_list.size(); ++i)
    if (class_list[i] != NULL)
      {
        class_list[i]->enter_my_subsection(prm);
        class_list[i]->declare_parameters(prm);
        class_list[i]->leave_my_subsection(prm);
      }
}


std::vector<std::string>
ParameterAcceptor::get_section_path() const
{
  Assert(acceptor_id < class_list.size(), ExcInternalError());
  std::vector<std::string> sections =
    Utilities::split_string_list(class_list[acceptor_id]->get_section_name(), sep);
  bool is_absolute = false;
  if (sections.size() > 1)
    {
      // Handle the cases of a leading "/"
      if (sections[0] == "")
        {
          is_absolute = true;
          sections.erase(sections.begin());
        }
    }
  if (is_absolute == false)
    {
      // In all other cases, we scan for earlier classes, and prepend the
      // first absolute path (in reverse order) we find to ours
      for (int i=acceptor_id-1; i>=0; --i)
        if (class_list[i] != NULL)
          if (class_list[i]->get_section_name().front() == sep)
            {
              bool has_trailing = class_list[i]->get_section_name().back() == sep;
              // Absolute path found
              auto secs = Utilities::split_string_list(class_list[i]->get_section_name(), sep);
              Assert(secs[0] == "", ExcInternalError());
              // Insert all sections except first and last
              sections.insert(sections.begin(), secs.begin()+1, secs.end()-(has_trailing ? 0 : 1));
              // exit from for cycle
              break;
            }
    }
  return sections;
}

void ParameterAcceptor::enter_my_subsection(ParameterHandler &prm)
{
  std::vector<std::string> sections = get_section_path();
  for (auto sec : sections)
    {
      prm.enter_subsection(sec);
    }
}

void ParameterAcceptor::leave_my_subsection(ParameterHandler &prm)
{
  std::vector<std::string> sections = get_section_path();
  for (auto sec : sections)
    {
      prm.leave_subsection();
    }
}



/// Conversion specializations.
/// string
template<>
std::shared_ptr<Patterns::PatternBase>  ParameterAcceptor::to_pattern<std::string>(const std::string &)
{
  return SP(new Patterns::Anything());
}

template<>
std::string ParameterAcceptor::to_string<std::string>(const std::string &entry)
{
  return entry;
}

template<>
std::string ParameterAcceptor::to_type<std::string>(const std::string &parameter)
{
  return parameter;
}

/// double
template<>
std::shared_ptr<Patterns::PatternBase>  ParameterAcceptor::to_pattern<double>(const double &)
{
  return SP(new Patterns::Double());
}

template<>
std::string ParameterAcceptor::to_string<double>(const double &entry)
{
  return std::to_string(entry);
}

template<>
double ParameterAcceptor::to_type<double>(const std::string &parameter)
{
  return std::stod(parameter);
}


/// int
template<>
std::shared_ptr<Patterns::PatternBase>  ParameterAcceptor::to_pattern<int>(const int &)
{
  return SP(new Patterns::Integer());
}

template<>
std::string ParameterAcceptor::to_string<int>(const int &entry)
{
  return std::to_string(entry);
}

template<>
int ParameterAcceptor::to_type<int>(const std::string &parameter)
{
  return std::stoi(parameter);
}


/// unsigned int
template<>
std::shared_ptr<Patterns::PatternBase>  ParameterAcceptor::to_pattern<unsigned int>(const unsigned int &)
{
  return SP(new Patterns::Integer(0));
}

template<>
std::string ParameterAcceptor::to_string<unsigned int>(const unsigned int &entry)
{
  return std::to_string(entry);
}

template<>
unsigned int ParameterAcceptor::to_type<unsigned int>(const std::string &parameter)
{
  return (unsigned int)std::stoi(parameter);
}


/// bool
template<>
std::shared_ptr<Patterns::PatternBase>  ParameterAcceptor::to_pattern<bool>(const bool &)
{
  return SP(new Patterns::Bool());
}

template<>
std::string ParameterAcceptor::to_string<bool>(const bool &entry)
{
  return std::string((entry ? "true" : "false"));
}

template<>
bool ParameterAcceptor::to_type<bool>(const std::string &parameter)
{
  return (parameter == "true" || parameter == "1");
}


/// point
#define MYP(dim) \
  template<>\
  std::shared_ptr<Patterns::PatternBase>  ParameterAcceptor::to_pattern<Point<dim> >(const Point<dim> &) {\
    return SP(new Patterns::List(Patterns::Double(),dim,dim));\
  }\
  \
  template<>\
  std::string ParameterAcceptor::to_string<Point<dim> >(const Point<dim> &entry) {\
    return print(entry);\
  }\
  \
  template<>\
  Point<dim> ParameterAcceptor::to_type<Point<dim> >(const std::string &parameter) {\
    Point<dim> p;\
    auto ps = Utilities::split_string_list(parameter);\
    AssertDimension(ps.size(), dim);\
    for(unsigned int i=0; i<dim; ++i)\
      p[i] = std::stod(ps[i]);\
    return p;\
  }

MYP(1)
MYP(2)
MYP(3)
#undef MYP


/// vector
#define MYV(type, sep) \
  template<>\
  std::shared_ptr<Patterns::PatternBase>  ParameterAcceptor::to_pattern<std::vector<type> >(const std::vector<type> &) {\
    type p;\
    return SP(new Patterns::List(*to_pattern(p),0,Patterns::List::max_int_value,sep));\
  }\
  \
  template<>\
  std::string ParameterAcceptor::to_string<std::vector<type> >(const std::vector<type> &entry) {\
    std::ostringstream s; \
    if(entry.size() > 0) \
      s << to_string(entry[0]); \
    for(unsigned int i=1; i<entry.size(); ++i) \
      s << std::string(sep) << " " << to_string(entry[i]); \
    return s.str();\
  }\
  \
  template<>\
  std::vector<type> ParameterAcceptor::to_type<std::vector<type> >(const std::string &parameter) {\
    auto ps = Utilities::split_string_list(parameter,sep[0]);\
    std::vector<type> p(ps.size());\
    for(unsigned int i=0; i<ps.size(); ++i)\
      p[i] = to_type<type>(ps[i]);\
    return p;\
  }

MYV(std::string, ",")
MYV(double, ",")
MYV(int, ",")
MYV(unsigned int, ",")
// MYV(bool, ",")
MYV(Point<1>, ";")
MYV(Point<2>, ";")
MYV(Point<3>, ";")


MYV(std::vector<std::string>, ";")
MYV(std::vector<double>, ";")
MYV(std::vector<int>, ";")
MYV(std::vector<unsigned int>, ";")
// MYV(std::vector<std::vector<bool> >, ";")
// MYV(std::vector<Point<1> >, ";")
// MYV(std::vector<Point<2> >, ";")
// MYV(std::vector<Point<3> >, ";")
#undef MYV

namespace
{
  template<class T>
  inline bool to_boost_any(const std::string &entry,
                           const Patterns::PatternBase &pattern,
                           boost::any &boost_parameter)
  {
    if (boost_parameter.type() == typeid(T *))
      {
        T &parameter = *(boost::any_cast<T *>(boost_parameter));
        parameter = ParameterAcceptor::to_type<T>(entry);
        AssertThrow(pattern.match(entry),
                    ExcMessage("The entry '"+entry+"' does not match"
                               " the pattern you specified: "+ pattern.description()));
        return true;
      }
    else
      {
        return false;
      }
  }
}

void ParameterAcceptor::parse_parameters(ParameterHandler &prm)
{
  for (auto &it : parameters)
    {
      const std::string &entry= prm.get(it.first);
      boost::any &boost_parameter = it.second;
      const Patterns::PatternBase &pattern = *patterns[it.first];

      if (to_boost_any<std::string>(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<double>(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<int>(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<unsigned int>(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<bool>(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<Point<1> >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<Point<2> >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<Point<3> >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<std::vector<std::string> >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<std::vector<double> >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<std::vector<int> >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<std::vector<unsigned int> >(entry, pattern, boost_parameter)) {}
//      else if (to_boost_any<std::vector<bool> >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<std::vector<Point<1> > >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<std::vector<Point<2> > >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<std::vector<Point<3> > >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<std::vector<std::vector<std::string> > >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<std::vector<std::vector<double> > >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<std::vector<std::vector<int> > >(entry, pattern, boost_parameter)) {}
      else if (to_boost_any<std::vector<std::vector<unsigned int> > >(entry, pattern, boost_parameter)) {}
//      else if (to_boost_any<std::vector<std::vector<bool> > >(entry, pattern, boost_parameter)) {}
      else
        {
          AssertThrow(false, ExcNotImplemented());
        }
    }
}

void ParameterAcceptor::parse_parameters_call_back() {}

D2K_NAMESPACE_CLOSE

