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
