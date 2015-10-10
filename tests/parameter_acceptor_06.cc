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
    add_parameter(prm, &my_int, "A parameter", "1", Patterns::Integer());
  };

  void log_info()
  {
    deallog << "My type: " << type(*this) << std::endl;
    deallog << my_int << std::endl;
  };

  unsigned int my_int;
};


int main ()
{
  initlog();
  Test<1> a;
  {
    // Since this class is constructed and immediately destructed,
    // it should not appear in the parameters below.
    Test<2> b;
    b.my_int = 0;
  }

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  a.log_info();

  ParameterAcceptor::log_info();
}
