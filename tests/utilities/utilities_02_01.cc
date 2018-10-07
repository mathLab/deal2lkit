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

#include <deal2lkit/utilities.h>

#include "../tests.h"


using namespace deal2lkit;


int
main()
{
  initlog();
  std::system("rm -rf ./tmp_dir*");
  deallog << "--> " << create_directory("tmp_dir000") << std::endl;
  std::string  new_dir = get_next_available_directory_name("tmp_dir", 3, 0, 10);
  unsigned int new_idx =
    get_next_available_index_directory_name("tmp_dir", 3, 0, 10);
  deallog << "--> " << dir_exists("tmp_dir001") << std::endl;
  deallog << "--> " << new_idx << std::endl;
  deallog << "--> " << new_dir << std::endl;
  deallog << "--> " << create_directory(new_dir) << std::endl;
  std::system("rm -rf ./tmp_dir*");
  deallog << "--> " << dir_exists("tmp_dir000") << std::endl;
  return 0;
}
