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

// Read a file in iges format, and write it out again in the same
// format.

#include "../tests.h"
#include <deal2lkit/parameter_acceptor.h>


using namespace deal2lkit;

template<int dim>
class Test : public deal2lkit::ParameterAcceptor
{
public:
  virtual void declare_parameters(ParameterHandler &prm)
  {
    prm.declare_entry("A double", "0.0", Patterns::Double(),
                      "Documentation");
  };

  virtual void parse_parameters(ParameterHandler &prm)
  {
    deallog << "Double: "
            << prm.get_double("A double") << std::endl;
  };
};


int main ()
{
  initlog();
  Test<2> a;
  Test<1> b;

  ParameterHandler prm;
  a.declare_all_parameters(prm);
  prm.log_parameters(deallog);
}
