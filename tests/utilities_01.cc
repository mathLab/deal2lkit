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

// Test if the "print" utilities function works as expected. This is
// somewhat the inverse of the dealii::Utilities::split_string_list(),
// but it works with arbitrary objects that support operator<<();

#include "tests.h"
#include "utilities.h"

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
