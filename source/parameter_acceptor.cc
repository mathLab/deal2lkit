#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/utilities.h>
#include <deal.II/base/point.h>
#include <deal.II/base/revision.h>
#include <fstream>

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
  prm.clear();
  declare_all_parameters(prm);
  if (filename != "")
    prm.read_input(filename);
  parse_all_parameters(prm);
  if (out_filename != "")
    {
      std::ofstream outfile(out_filename.c_str());
      Assert(outfile, ExcIO());
      outfile << "# Parameter file generated with " << std::endl
              << "# D2K_GIT_BRANCH=       " << D2K_GIT_BRANCH << std::endl
              << "# D2K_GIT_SHORTREV=     " << D2K_GIT_SHORTREV << std::endl
              << "# DEAL_II_GIT_BRANCH=   " << DEAL_II_GIT_BRANCH  << std::endl
              << "# DEAL_II_GIT_SHORTREV= " << DEAL_II_GIT_SHORTREV << std::endl;
      prm.print_parameters(outfile, ParameterHandler::ShortText);
    }
}

void
ParameterAcceptor::clear()
{
  class_list.clear();
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
        prm.enter_subsection(class_list[i]->get_section_name());
        class_list[i]->parse_parameters(prm);
        class_list[i]->parse_parameters_call_back();
        prm.leave_subsection();
      }
}

void ParameterAcceptor::declare_all_parameters(ParameterHandler &prm)
{
  for (unsigned int i=0; i< class_list.size(); ++i)
    if (class_list[i] != NULL)
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
      else if (it->second.type() == typeid(std::vector<unsigned int> *))
        {
          std::vector<unsigned int> &int_list = *(boost::any_cast<std::vector<unsigned int>*>(it->second));
          std::vector<std::string> string_list = Utilities::split_string_list(prm.get(it->first));
          int_list.resize(string_list.size());
          for (unsigned int i=0; i<string_list.size(); ++i)
            {
              std::istringstream reader(string_list[i]);
              reader >> int_list[i];// = std::stoul(string_list[i]);
            }
        }
      else if (it->second.type() == typeid(std::vector<double> *))
        {
          std::vector<double> &double_list = *(boost::any_cast<std::vector<double>*>(it->second));
          double_list = Utilities::string_to_double(Utilities::split_string_list(prm.get(it->first)));
        }
      else if (it->second.type() == typeid(std::vector<std::vector<unsigned int> > *))
        {
          std::vector<std::vector<unsigned int> > &int_table =
            *(boost::any_cast<std::vector<std::vector<unsigned int> >*>(it->second));

          std::vector<std::string> string_list_of_lists = Utilities::split_string_list(prm.get(it->first), ';');
          int_table.resize(string_list_of_lists.size());
          for (unsigned int i=0; i<int_table.size(); ++i)
            {
              std::vector<std::string> string_list = Utilities::split_string_list(string_list_of_lists[i], ',');
              int_table[i].resize(string_list.size());
              for (unsigned int j=0; j<int_table[i].size(); ++j)
                {
                  std::istringstream reader(string_list[j]);
                  reader >> int_table[i][j];
                }
            }
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
        }
    }
}


void ParameterAcceptor::parse_parameters_call_back() {}
