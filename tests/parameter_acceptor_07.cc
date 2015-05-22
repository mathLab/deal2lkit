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

// Test parmeter table

#include "tests.h"
#include "utilities.h"
#include "parameter_acceptor.h"

template<int dim>
class Test : public ParameterAcceptor
{
public:
  virtual void declare_parameters(ParameterHandler &prm)
  {
    add_parameter(prm, &my_list, "A table", "1,2,3;5,6 ; 7,8,9, 10", Patterns::Anything());
  };

  void log_info()
  {
    deallog << "My type: " << type(*this) << std::endl;
    for (unsigned int i=0; i<my_list.size(); ++i)
      deallog << i << ": " << print(my_list[i]) << std::endl;
  };

  std::vector<std::vector<unsigned int> > my_list;
};


int main ()
{
  initlog();
  Test<1> a;

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);
  a.log_info();

}
