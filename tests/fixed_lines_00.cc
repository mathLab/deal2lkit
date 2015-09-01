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

#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "tests.h"
#include "utilities.h"

int main ()
{
  initlog();

  TimeUtilities tu;

  unsigned int  N = 3;
  int           a = 100;
  double        b = 3.14;

  OverWriteStream<dealii::LogStream> out(N, deallog);

  out << a;
  tu.sleep(500);
  out << std::endl;              // manipulator inserted alone
  tu.sleep(500);
  out << b << ", this is a longer line" << std::endl;
  tu.sleep(500);
  out << a * b << std::endl;  // manipulator in concatenated insertion
  tu.sleep(500);

  out.clear(true);

  for (unsigned int i=0; i<100; ++i)
    {
      out << "Line " << std::setw(6) << i  << std::endl;
      tu.sleep(100);
    }

  return 0;
}
