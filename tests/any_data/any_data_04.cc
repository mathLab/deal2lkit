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

#include "../tests.h"
#include <deal2lkit/any_data.h>
#include <deal.II/base/tensor.h>



using namespace deal2lkit;


int main ()
{
  initlog();

  AnyData data;

  data.add_copy(0, "zero");
  deallog << (data.have("zero") ? "OK" : "KO") << std::endl;

}
