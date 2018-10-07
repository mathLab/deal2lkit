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

#include <stdlib.h>

#include "../tests.h"



int main()
{
  initlog();
  const std::string name_dir = "directory_test";
  std::string       cmd      = "";

  cmd = name_dir;
  if (std::system(cmd.c_str()))
    {
      cmd = "test -d " + name_dir;
      if (std::system(cmd.c_str()))
        deallog << "The directory " + name_dir + " has been correctly created."
                << std::endl;

      cmd = "rmdir" + name_dir;
      if (std::system(cmd.c_str()))
        {
          deallog << "The directory " + name_dir +
                       " has been correctly removed."
                  << std::endl;
        }
      else
        {
          deallog << "Impossible to remove " + name_dir + "." << std::endl;
        }
    }
  else
    {
      deallog << "Impossible to create " + name_dir + "." << std::endl;
    }
}
