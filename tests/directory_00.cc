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
#include <stdlib.h>

int main ()
{
  initlog();
  const std::string name_dir = "directory_test";
  std::string cmd = "";

  cmd = name_dir;
  if (std::system(cmd.c_str()))
    {
      cmd = "test -d "+name_dir;
      if (std::system(cmd.c_str()))
        deallog << "The directory " + name_dir + " has been correctly created." << std::endl;

      cmd = "rmdir" + name_dir;
      if (std::system(cmd.c_str()))
        {
          deallog << "The directory " + name_dir + " has been correctly removed." << std::endl;
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
