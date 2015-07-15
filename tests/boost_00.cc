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
#ifdef DEAL_II_SAK_WITH_BOOST
  deallog << "Boost library found " << std::endl;
#else
  deallog << "Boost library not found " << std::endl;
#endif
  deallog << "Program works " << std::endl;
}
