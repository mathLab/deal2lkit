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
#include "parameter_acceptor.h"
#include "parsed_mapped_functions.h"


int main ()
{
  initlog();
  ParsedMappedFunctions<2,3> dirichlet("Dirichlet BCs", "u,u,p", "0:0;1,1:2,6:0;1;2","0:x;y;0, 1:0;0;0,6:y*k;0;k","k=1");

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  Point<2> p(2,3);

  deallog << "Component mask id 0: " <<  dirichlet.get_mapped_mask(0) << std::endl;
  deallog << "Component mask id 1: " <<  dirichlet.get_mapped_mask(1) << std::endl;
  deallog << "Component mask id 6: " <<  dirichlet.get_mapped_mask(6) << std::endl;
  deallog << "Parsed Function on id 0: " <<  (*dirichlet.get_mapped_function(0)).value(p) << std::endl;
  deallog << "Parsed Function on id 1: " <<  dirichlet.get_mapped_function(1)->value(p) << std::endl;
  deallog << "Parsed Function on id 6: " <<  dirichlet.get_mapped_function(6)->value(p) << std::endl;
}
