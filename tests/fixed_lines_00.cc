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
#include <stdlib.h>
#include <iostream>

int main ()
{
  initlog();

  fixed_lines stream_out_test(3, std::cout);

  std::string testo1 = "testo1";
  std::string testo2 = "testo2";

  stream_out_test.print_line(testo1);
  stream_out_test.print_line(testo1);
  stream_out_test.print_line(testo1);
  stream_out_test.print_line(testo2);

  // std::cout << "\n\n\n\n";
}
