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
#include "sak_data.h"
#include <deal.II/base/tensor.h>


template<int dim, int spacedim>
class Energy
{
public:
  template <typename Number>
  void fill_data(SAKData &d)
  {
    Tensor<1,spacedim, Number> p;
    p[0] = 1.0;
    std::string suffix = demangle(typeid(Number).name());
    d.add_copy(p, "u"+suffix);
  }

  template <typename Number>
  Number energy(const SAKData &d) const
  {
    std::string suffix = demangle(typeid(Number).name());
    auto &p = d.get<Tensor<1,spacedim,Number> >("u"+suffix);
    return p.norm();
  }
};

template<int dim, int spacedim, class EnergyClass>
class Problem
{
public:
  EnergyClass e;

  void run()
  {
    SAKData d;
    e.template fill_data<double>(d);
    e.template fill_data<int>(d);
    deallog << "Double: " <<  e.template energy<double>(d) << std::endl;
    deallog << "Int: " << e.template energy<int>(d) << std::endl;
    d.print_info(deallog);
  }
};


int main ()
{
  initlog();

  Problem<1,1,Energy<1,1> > pb;
  pb.run();
}
