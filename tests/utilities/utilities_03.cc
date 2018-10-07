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

// Test if the "copy files" utilities function works as expected

#include <deal2lkit/utilities.h>

#include "../tests.h"


using namespace deal2lkit;

int main()
{
  initlog();
  std::system("touch test1.txt");
  std::system("touch test2.txt");
  std::system("rm -rf new_folder*");
  deallog << "copy   --> " << copy_file("test1.txt", "tmp") << std::endl;
  deallog << "exist  --> " << file_exists("test1.txt") << std::endl;
  deallog << "exist  --> " << file_exists("test2.txt") << std::endl;
  std::string folder = get_next_available_directory_name("new_folder", 3);
  deallog << "create --> " << folder << " - " << create_directory(folder)
          << std::endl;
  deallog << "copy   --> " << copy_files("test1.txt test2.txt", folder)
          << std::endl;
  deallog << "exist  --> " << file_exists(folder + "/" + "test1.txt")
          << std::endl;
  deallog << "exist  --> " << file_exists(folder + "/" + "test2.txt")
          << std::endl;
  std::system("rm -rf new_folder*");
  deallog << "exist  --> " << dir_exists("new_folder000") << std::endl;
}
