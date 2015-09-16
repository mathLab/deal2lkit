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

#include <deal.II/base/point.h>

template<int dim>
class Test : public ParameterAcceptor
{
public:
  virtual void declare_parameters(ParameterHandler &prm)
  {
    std::string def = "0.";
    for (int i=1; i<dim; ++i)
      def += ",0.";
    add_parameter(prm, &p, "A point", def, Patterns::List(Patterns::Double(), dim, dim));
  };

  void log_info()
  {
    deallog << "My type: " << type(*this) << std::endl
            << "p: " << p << std::endl;
  }

private:
  Point<dim> p;
};


int main ()
{
  initlog();
  Test<1> a;
  Test<2> b;
  Test<3> c;

  ParameterHandler prm;
  ParameterAcceptor::declare_all_parameters(prm);
  prm.read_input_from_string(""
                             "subsection Test<1>\n"
                             "  set A point = 1.0\n"
                             "end\n"
                             "subsection Test<2>\n"
                             "  set A point = 1.0, 2.0\n"
                             "end\n"
                             "subsection Test<3>\n"
                             "  set A point = 1.0, 2.0, 3.0\n"
                             "end\n");

  prm.log_parameters(deallog);
  ParameterAcceptor::parse_all_parameters(prm);

  a.log_info();
  b.log_info();
  c.log_info();
}
