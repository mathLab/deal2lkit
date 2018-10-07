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



#include <deal.II/base/path_search.h>

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/utilities.h>

#include "../tests.h"


using namespace deal2lkit;

class MyClass : public ParameterAcceptor
{
public:
  MyClass(const std::string &name = "My Class", const unsigned int &i = 0) :
    ParameterAcceptor(name),
    i(i)
  {
    add_parameter(this->i, "An integer");
  };

private:
  unsigned int i;
};

int main()
{
  initlog();

  MyClass f0("One", 1);
  MyClass f1("Two", 2);
  MyClass f2("Three", 3);

  ParameterAcceptor::prm.log_parameters(deallog);
}
