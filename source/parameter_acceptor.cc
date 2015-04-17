#include "parameter_acceptor.h"
#include "utilities.h"

#include <deal.II/base/point.h>

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
ParameterAcceptor::initialize(const std::string filename)
{
  prm.clear();
  declare_all_parameters(prm);
  if (filename != "")
    prm.read_input(filename);
  parse_all_parameters(prm);
}

void
ParameterAcceptor::clear()
{
  class_list.clear();
}

void ParameterAcceptor::parse_all_parameters(ParameterHandler &prm)
{
  for (unsigned int i=0; i< class_list.size(); ++i)
    {
      prm.enter_subsection(class_list[i]->get_section_name());
      class_list[i]->parse_parameters(prm);
      prm.leave_subsection();
    }
}

void ParameterAcceptor::declare_all_parameters(ParameterHandler &prm)
{
  for (unsigned int i=0; i< class_list.size(); ++i)
    {
      prm.enter_subsection(class_list[i]->get_section_name());
      class_list[i]->declare_parameters(prm);
      prm.leave_subsection();
    }
}

void ParameterAcceptor::parse_parameters(ParameterHandler &prm)
{
  for (auto it = parameters.begin(); it != parameters.end(); ++it)
    {
      if (it->second.type() == typeid(std::string *))
        {
          *(boost::any_cast<std::string *>(it->second)) = prm.get(it->first);
        }
      else if (it->second.type() == typeid(double *))
        {
          *(boost::any_cast<double *>(it->second)) = prm.get_double(it->first);
        }
      else if (it->second.type() == typeid(int *))
        {
          *(boost::any_cast<int *>(it->second)) = prm.get_integer(it->first);
        }
      else if (it->second.type() == typeid(unsigned int *))
        {
          *(boost::any_cast<unsigned int *>(it->second)) = prm.get_integer(it->first);
        }
      else if (it->second.type() == typeid(bool *))
        {
          *(boost::any_cast<bool *>(it->second)) = prm.get_bool(it->first);
        }
      // Here we have the difficult types...
      // First all point types.
      else if (it->second.type() == typeid(Point<1> *))
        {
          std::vector<double> p =
            Utilities::string_to_double(Utilities::split_string_list(prm.get(it->first)));
          AssertDimension(p.size(), 1);

          Point<1> &pp = *(boost::any_cast<Point<1>*>(it->second));
          pp[0] = p[0];
        }
      else if (it->second.type() == typeid(Point<2> *))
        {
          std::vector<double> p =
            Utilities::string_to_double(Utilities::split_string_list(prm.get(it->first)));
          AssertDimension(p.size(), 2);

          Point<2> &pp = *(boost::any_cast<Point<2>*>(it->second));
          pp[0] = p[0];
          pp[1] = p[1];
        }
      else if (it->second.type() == typeid(Point<3> *))
        {
          std::vector<double> p =
            Utilities::string_to_double(Utilities::split_string_list(prm.get(it->first)));
          AssertDimension(p.size(), 3);

          Point<3> &pp = *(boost::any_cast<Point<3>*>(it->second));
          pp[0] = p[0];
          pp[1] = p[1];
          pp[2] = p[2];
        }
      else if (it->second.type() == typeid(std::vector<std::string> *))
        {
          std::vector<std::string> &string_list = *(boost::any_cast<std::vector<std::string>*>(it->second));
          string_list = Utilities::split_string_list(prm.get(it->first));
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
        }
    }
}


void ParameterAcceptor::parse_parameters_call_back() {}
