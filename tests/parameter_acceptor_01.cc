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

// Read a file in iges format, and write it out again in the same
// format.

#include "tests.h"
#include "parameter_acceptor.h"

template<int dim>
class Test : public ParameterAcceptor
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
