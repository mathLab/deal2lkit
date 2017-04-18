//-----------------------------------------------------------
//
//    Copyright (C) 2016 by the deal2lkit authors
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

// test the to_double() function


#include "../tests.h"

#include <fstream>

#include <deal2lkit/sacado_tools.h>


using namespace deal2lkit;

typedef Sacado::Fad::DFad<double>  Sdouble;
typedef Sacado::Fad::DFad<Sdouble> SSdouble;


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  deallog << " TEST : " << std::endl;
  double num = 4;

  deallog << " double " << std::endl;
  deallog << "  Square root of " << num
          << " is " << std::sqrt(num) << std::endl;
  deallog << "  sqrt(x)^2  = "
          << std::sqrt(num)*std::sqrt(num)
          << std::endl;

  deallog << " Sdouble " << std::endl;
  Sdouble Snum(2,1,4.0);
  deallog << "  Square root of " << Snum
          << " is " << std::sqrt(Snum) << std::endl;
  deallog << "  sqrt(x)^2  = "
          << std::sqrt(Snum)*std::sqrt(Snum)
          << std::endl;

  deallog << " SSdouble " << std::endl;
  SSdouble SSnum(2,1,Snum);
  deallog << "  Square root of " << SSnum
          << " is " << std::sqrt(SSnum) << std::endl;
  deallog << "  sqrt(x)^2  = "
          << std::sqrt(SSnum)*std::sqrt(SSnum)
          << std::endl;
}
