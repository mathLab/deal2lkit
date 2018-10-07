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


#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_finite_element.h>
#include <deal2lkit/utilities.h>

#include "../tests.h"


using namespace deal2lkit;


int main()
{
  initlog();
  ParsedFiniteElement<2, 2> fe_builder(
    "FE", "FESystem[FE_Q(2)^d-FE_DGP(1)]", "u,u,p");

  ParameterAcceptor::initialize();
  auto fe = fe_builder();

  deallog << "N components: " << fe_builder.n_components()
          << ", N blocks: " << fe_builder.n_blocks() << std::endl
          << "Components: " << fe_builder.get_component_names()
          << ", Blocks: " << fe_builder.get_block_names() << std::endl
          << " position of u: " << fe_builder.get_first_occurence("u")
          << std::endl
          << " position of p: " << fe_builder.get_first_occurence("p")
          << std::endl
          << " is u vector?: " << fe_builder.is_vector("u") << std::endl
          << " is p vector?: " << fe_builder.is_vector("p") << std::endl;
  std::vector<unsigned int> b = fe_builder.get_component_blocks();
  deallog << "Component blocks: " << b[0];
  for (unsigned int i = 1; i < fe_builder.n_components(); ++i)
    deallog << ", " << b[i];
  deallog << std::endl;
}
