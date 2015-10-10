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

#include "tests.h"
#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>


using namespace deal2lkit;

template<int dim>
class Test : public ParameterAcceptor
{
public:
  virtual void declare_parameters(ParameterHandler &prm)
  {
    add_parameter(prm, &string_list, "A list of strings", "ciao, mamma", Patterns::List(Patterns::Anything()));
  };

  void log_info()
  {
    deallog << "My type: " << type(*this) << std::endl;
    for (unsigned int i=0; i<string_list.size(); ++i)
      deallog << string_list[i] << std::endl;
  };

private:
  std::vector<std::string> string_list;
};


int main ()
{
  initlog();
  Test<1> a;

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  a.log_info();
}
