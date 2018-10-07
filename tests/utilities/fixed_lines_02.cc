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

#include <deal2lkit/utilities.h>
#include <stdlib.h>

#include <iomanip>
#include <iostream>

#include "../tests.h"


using namespace deal2lkit;

int
main()
{
  std::fstream fs("test_fstream.txt", std::fstream::out);

  TimeUtilities tu;

  unsigned int N = 3;
  int          a = 100;
  double       b = 3.14;

  OverWriteStream<> out(N, fs);

  out << a;
  tu.sleep(500);
  out << std::endl; // manipulator inserted alone
  tu.sleep(500);
  out << b << ", this is a longer line" << std::endl;
  tu.sleep(500);
  out << a * b << std::endl; // manipulator in concatenated insertion
  tu.sleep(500);

  out.clear(true);

  for (unsigned int i = 0; i < 100; ++i)
    {
      out << "Line " << std::setw(6) << i << std::endl;
      tu.sleep(100);
    }

  out.end();
  fs.close();
  rename_file("test_fstream.txt", "output");
  return 0;
}
