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
#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>

template<int dim>
class Test : public ParameterAcceptor
{
public:
  virtual void declare_parameters(ParameterHandler &prm)
  {
    add_parameter(prm, &int_list, "A list of ints", "1, 2", Patterns::List(Patterns::Integer()));
  };

  void log_info()
  {
    deallog << "My type: " << type(*this) << std::endl;
    for (unsigned int i=0; i<int_list.size(); ++i)
      deallog << int_list[i] << std::endl;
  };

private:
  std::vector<unsigned int> int_list;
};


int main ()
{
  initlog();
  Test<1> a;

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  a.log_info();
}
