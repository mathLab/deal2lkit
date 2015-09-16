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
