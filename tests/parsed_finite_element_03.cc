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
#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_finite_element.h>


int main ()
{
  initlog();
  ParsedFiniteElement<2,2> fe_builder("FE", "FESystem[FE_Q(2)^d-FE_DGP(1)]",
                                      "u,u,p");

  ParameterAcceptor::initialize();
  auto fe = fe_builder();

  deallog << "N components: " << fe_builder.n_components()
          << ", N blocks: " << fe_builder.n_blocks() << std::endl
          << "Components: " << fe_builder.get_component_names()
          << ", Blocks: " << fe_builder.get_block_names() << std::endl;
  std::vector<unsigned int> b = fe_builder.get_component_blocks();
  deallog << "Component blocks: " << b[0];
  for (unsigned int i=1; i<fe_builder.n_components(); ++i)
    deallog << ", " << b[i];
  deallog << std::endl;
}
