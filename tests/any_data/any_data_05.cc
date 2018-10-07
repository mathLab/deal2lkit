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

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal2lkit/any_data.h>
#include <deal2lkit/utilities.h>

#include "../tests.h"


using namespace deal2lkit;


template <int dim, int spacedim>
class Energy
{
public:
  template <typename Number>
  void fill_data(AnyData &d)
  {
    Tensor<1, spacedim, Number> p;
    p[0]               = 1.0;
    std::string suffix = demangle(typeid(Number).name());
    d.add_copy(p, "u" + suffix);
  }

  template <typename Number>
  Number energy(const AnyData &d) const
  {
    std::string suffix = demangle(typeid(Number).name());
    auto &      p      = d.get<Tensor<1, spacedim, Number>>("u" + suffix);
    return p.norm();
  }
};

template <int dim, int spacedim, class EnergyClass>
class Problem
{
public:
  EnergyClass e;

  void run()
  {
    AnyData d;
    e.template fill_data<double>(d);
    e.template fill_data<int>(d);
    deallog << "Double: " << e.template energy<double>(d) << std::endl;
    deallog << "Int: " << e.template energy<int>(d) << std::endl;
    d.print_info(deallog);
  }
};


int main()
{
  initlog();

  Problem<1, 1, Energy<1, 1>> pb;
  pb.run();
}
