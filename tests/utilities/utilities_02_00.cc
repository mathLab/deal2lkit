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


// Test if the "create directory" utilities function works as expected.

#include "../tests.h"
#include <deal2lkit/utilities.h>



using namespace deal2lkit;

int main ()
{
  initlog();

  std::system("rm -rf dsfas*");
  deallog << "--> " << create_directory("dsfas000") << std::endl;
  deallog << "--> " << create_directory("dsfas001") << std::endl;
  deallog << "--> " << create_directory("dsfas002") << std::endl;
  deallog << "--> " << dir_exists("dsfas002") << std::endl;
  std::system("rm -rf dsfas*");
  deallog << "--> " << dir_exists("dsfas002") << std::endl;
  return 0;
}
