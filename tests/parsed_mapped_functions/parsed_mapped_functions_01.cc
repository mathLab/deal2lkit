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

// test basic functionalities


#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_mapped_functions.h>
#include <deal2lkit/utilities.h>

#include "../tests.h"


using namespace deal2lkit;


int
main()
{
  initlog();
  ParsedMappedFunctions<2> pmf("Mapped Functions",
                               3,
                               "u,u,p",
                               "0=0;1 % 1=2 % 6=0;1;2",
                               "0=x;y;0 % 1=0;0;0 % 6=y*k;beta*y;k",
                               "k=1,beta=2");

  dealii::ParameterAcceptor::initialize(
    SOURCE_DIR "/parameters/parsed_mapped_functions_01.prm",
    "used_parameters.prm");
  dealii::ParameterAcceptor::prm.log_parameters(deallog);

  Point<2> p(2, 3);

  std::vector<unsigned int> ids(3);
  for (unsigned int i = 0; i < ids.size(); ++i)
    ids[i] = i;

  unsigned int id = 3;

  if (std::find(ids.begin(), ids.end(), id) != ids.end())
    deallog << "c'e" << std::endl;


  deallog << "Component mask id 0: " << pmf.get_mapped_mask(0) << std::endl;
  deallog << "Component mask id 1: " << pmf.get_mapped_mask(1) << std::endl;
  deallog << "Component mask id 6: " << pmf.get_mapped_mask(6) << std::endl;
  deallog << "Parsed Function on id 0: "
          << (*pmf.get_mapped_function(0)).value(p) << std::endl;
  deallog << "Parsed Function on id 1: " << pmf.get_mapped_function(1)->value(p)
          << std::endl;
  deallog << "Parsed Function on id 6: " << pmf.get_mapped_function(6)->value(p)
          << std::endl;
}
