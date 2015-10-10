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

#include "tests.h"
#include <deal2lkit/utilities.h>


using namespace deal2lkit;

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
