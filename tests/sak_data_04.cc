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

  data.add_copy(0, "zero");
  deallog << (data.have("zero") ? "OK" : "KO") << std::endl;

}
