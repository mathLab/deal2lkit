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

#include <deal2lkit/any_data.h>

#include "../tests.h"



using namespace deal2lkit;


int main()
{
  initlog();

  AnyData data;


  const unsigned int n_q = 5;

  std::vector<double> v_double(n_q);
  std::vector<int>    v_int(n_q);

  data.add_ref(v_int, "int_ref");
  data.add_copy(v_double, "double_copy");

  AnyData newdata = data;

  deallog << "verify that add_ref actually stored the reference" << std::endl;

  v_int[0]    = 7.0;
  v_double[0] = 7.0;

  deallog << "double_copy: ";
  for (unsigned int i = 0; i < n_q; ++i)
    {
      deallog << data.get<std::vector<double>>("double_copy")[i];
    }
  deallog << std::endl;
  deallog << "newdata::double_copy: ";
  for (unsigned int i = 0; i < n_q; ++i)
    {
      deallog << newdata.get<std::vector<double>>("double_copy")[i];
    }
  deallog << std::endl;
  deallog << "int_ref: ";
  for (unsigned int i = 0; i < n_q; ++i)
    deallog << data.get<std::vector<int>>("int_ref")[i] << " ";
  deallog << std::endl;
  deallog << "newdata::int_ref: ";
  for (unsigned int i = 0; i < n_q; ++i)
    deallog << newdata.get<std::vector<int>>("int_ref")[i] << " ";
  deallog << std::endl;
  deallog << std::endl;

  // How to manage the stored data
  std::vector<int> &   vi = data.get<std::vector<int>>("int_ref");
  std::vector<double> &vd = data.get<std::vector<double>>("double_copy");
  for (unsigned int i = 0; i < n_q; ++i)
    {
      vi[i] += i * i;
      vd[i] += i * i;
    }

  deallog << "read the stored data" << std::endl;
  deallog << "double: ";
  for (unsigned int i = 0; i < n_q; ++i)
    deallog << data.get<std::vector<double>>("double_copy")[i] << " ";
  deallog << std::endl;
  deallog << "int: ";
  for (unsigned int i = 0; i < n_q; ++i)
    deallog << data.get<std::vector<int>>("int_ref")[i] << " ";
  deallog << std::endl;

  deallog << "read newdata" << std::endl;
  deallog << "double: ";
  for (unsigned int i = 0; i < n_q; ++i)
    deallog << newdata.get<std::vector<double>>("double_copy")[i] << " ";
  deallog << std::endl;
  deallog << "int: ";
  for (unsigned int i = 0; i < n_q; ++i)
    deallog << newdata.get<std::vector<int>>("int_ref")[i] << " ";
  deallog << std::endl;
}
