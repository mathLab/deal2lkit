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

// test acts_on_this_id


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
                               "0=u % 1=1 % 6=u;p",
                               "0=x;y;0 % 1=0;0;0 % 6=y*k;0;k",
                               "k=1");

  dealii::ParameterAcceptor::initialize();
  dealii::ParameterAcceptor::prm.log_parameters(deallog);

  unsigned int id_true  = 0;
  unsigned int id_false = 7;

  deallog << "Mapped function acts on id 0 (expected true) --> "
          << pmf.acts_on_id(id_true) << std::endl;
  deallog << "Mapped function acts on id 7 (expected false) --> "
          << pmf.acts_on_id(id_false) << std::endl;
}
