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

// Test if the "print" utilities function works as expected. This is
// somewhat the inverse of the dealii::Utilities::split_string_list(),
// but it works with arbitrary objects that support operator<<();

#include "../tests.h"
#include <deal2lkit/utilities.h>


using namespace deal2lkit;

int main ()
{
  initlog();

  std::vector<unsigned int> t(4, 1);
  t[2] = 3;

  deallog << print(t) << std::endl;
  deallog << print(t, ";") << std::endl;

  std::vector<std::string> s(4, "u");
  s[3] = "p";

  deallog << print(s) << std::endl;
  deallog << print(s, "+") << std::endl;

}
