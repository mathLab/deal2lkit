//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------


#include "tests.h"
#include "utilities.h"
#include "parameter_acceptor.h"

template<int dim>
class Test : public ParameterAcceptor
{
public:
  virtual void declare_parameters(ParameterHandler &prm)
  {
    add_parameter(prm, &a, "A double", "1.0", Patterns::Double());
    add_parameter(prm, &b, "An int", "2", Patterns::Integer());
    add_parameter(prm, &c, "A string", "Ciao", Patterns::Anything());
    add_parameter(prm, &d, "A bool", "true", Patterns::Bool());
  };

  void log_info()
  {
    deallog << "My type: " << type(*this) << std::endl
            << "a: " << a << std::endl
            << "b: " << b << std::endl
            << "c: " << c << std::endl
            << "d: " << (d ? "true" : "false") << std::endl;
  }

private:
  double a;
  int b;
  std::string c;
  bool d;

};


int main ()
{
  initlog();
  Test<1> a;
  Test<2> b;

  ParameterHandler prm;
  ParameterAcceptor::declare_all_parameters(prm);
  ParameterAcceptor::parse_all_parameters(prm);
  prm.log_parameters(deallog);

  a.log_info();
  b.log_info();

  prm.read_input_from_string(""
                             "subsection Test<1>\n"
                             "  set A double = 3.0\n"
                             "end\n"
                             "subsection Test<2>\n"
                             "  set A double = 5.0\n"
                             "end\n");

  prm.log_parameters(deallog);
  ParameterAcceptor::parse_all_parameters(prm);

  a.log_info();
  b.log_info();

}
