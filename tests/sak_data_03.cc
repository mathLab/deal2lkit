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



int main ()
{
  initlog();

  SAKData data;


  const unsigned int n_q = 5;

  std::vector<double> v_double(n_q);
  std::vector<int> v_int(n_q);

  data.add_ref<std::vector<int> >(v_int, "int_ref");
  data.add_copy<std::vector<double> >(v_double, "double_copy");

  SAKData newdata = data;

  deallog << "verify that add_ref actually stored the reference"
          << std::endl;

  v_int[0] = 7.0;
  v_double[0] = 7.0;

  deallog << "double_copy: ";
  for (unsigned int i=0; i<n_q; ++i)
    {
      deallog << data.get<std::vector<double> > ("double_copy")[i];
    }
  deallog << std::endl;
  deallog << "newdata::double_copy: ";
  for (unsigned int i=0; i<n_q; ++i)
    {
      deallog << newdata.get<std::vector<double> > ("double_copy")[i];
    }
  deallog << std::endl;
  deallog << "int_ref: ";
  for (unsigned int i=0; i<n_q; ++i)
    deallog << data.get<std::vector<int> >("int_ref")[i] << " ";
  deallog << std::endl;
  deallog << "newdata::int_ref: ";
  for (unsigned int i=0; i<n_q; ++i)
    deallog << newdata.get<std::vector<int> >("int_ref")[i] << " ";
  deallog << std::endl;
  deallog << std::endl;

//How to manage the stored data
  auto &vi = data.get<std::vector<int> >("int_ref");
  auto &vd = data.get<std::vector<double> >("double_copy");
  for (unsigned int i=0; i<n_q; ++i)
    {
      vi[i] += i*i;
      vd[i] += i*i;
    }

  deallog << "read the stored data" << std::endl;
  deallog << "double: ";
  for (unsigned int i=0; i<n_q; ++i)
    deallog << data.get<std::vector<double> >("double_copy")[i] << " ";
  deallog << std::endl;
  deallog << "int: ";
  for (unsigned int i=0; i<n_q; ++i)
    deallog << data.get<std::vector<int> >("int_ref")[i] << " ";
  deallog << std::endl;

  deallog << "read newdata" << std::endl;
  deallog << "double: ";
  for (unsigned int i=0; i<n_q; ++i)
    deallog << newdata.get<std::vector<double> >("double_copy")[i] << " ";
  deallog << std::endl;
  deallog << "int: ";
  for (unsigned int i=0; i<n_q; ++i)
    deallog << newdata.get<std::vector<int> >("int_ref")[i] << " ";
  deallog << std::endl;
}
