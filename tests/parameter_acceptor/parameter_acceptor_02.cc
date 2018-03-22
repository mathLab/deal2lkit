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

#include "../tests.h"
#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>


using namespace deal2lkit;

template<int dim>
class Test : public deal2lkit::ParameterAcceptor
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
  deal2lkit::ParameterAcceptor::declare_all_parameters(prm);
  deal2lkit::ParameterAcceptor::parse_all_parameters(prm);
  prm.log_parameters(deallog);

  a.log_info();
  b.log_info();

  prm.parse_input_from_string(""
                              "subsection Test<1>\n"
                              "  set A double = 3.0\n"
                              "end\n"
                              "subsection Test<2>\n"
                              "  set A double = 5.0\n"
                              "end\n");

  prm.log_parameters(deallog);
  deal2lkit::ParameterAcceptor::parse_all_parameters(prm);

  a.log_info();
  b.log_info();

}
