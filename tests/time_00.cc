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

// This test is for the TimeUtilities class of utilities.h

#include "tests.h"
#include <deal2lkit/utilities.h>
#include <chrono>
#include <thread>

int main ()
{
  initlog();
  clock_t t;

  TimeUtilities time_utilities;

  deallog << " START " << std::endl;
  for (unsigned int i = 2; i < 6; ++i)
    {
      time_utilities.get_start_time();
      time_utilities.sleep(i*1000);
      time_utilities.get_end_time();

      deallog << " STEP:  " << time_utilities.get_num_measures() << std::endl;
      deallog << " SLEEP: " << int(i)                    << " sec" << std::endl;
      deallog << " TIMER: " << int(time_utilities[i-2])  << " sec" << std::endl;

      if ( abs(time_utilities[i-2] - i) < 1)
        {
          deallog << " OK " << std::endl;
        }
      else
        {
          deallog << " FAIL " << std::endl;
        }

    }
  deallog << " END " << std::endl;
}
