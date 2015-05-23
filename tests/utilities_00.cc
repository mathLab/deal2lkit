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
#include "utilities.h"

int main ()
{
  initlog();

  std::vector<std::string> components;
  components.push_back("u");
  components.push_back("u");
  components.push_back("p");

  std::vector<std::string>::iterator it;

  deallog << "original vector contains:" << std::endl;
  for (it = components.begin(); it!=components.end(); ++it)
    deallog << *it << std::endl;

  std::vector<std::string> uvec = unique(components);

  deallog << "unique vector contains:" << std::endl;
  for (it = uvec.begin(); it !=uvec.end(); ++it)
    deallog << *it << std::endl;
}
