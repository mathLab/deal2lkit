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

#include <stdio.h>
// #include <sys/types.h> // for ftruncate...



int
main ()
{
  fputs("output1\n",stdout);
  fputs("output2\n",stdout);
  fputs("\033[A\033[2K\033[A\033[2K",stdout);
  rewind(stdout);
  // ftruncate(1,0);
  fputs("output3\n",stdout);
  fputs("output4\n",stdout);
  return 0;
}
