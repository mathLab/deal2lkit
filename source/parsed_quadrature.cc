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

#include <deal2lkit/parsed_quadrature.h>

D2K_NAMESPACE_OPEN

template <int dim>
ParsedQuadrature<dim>::
ParsedQuadrature(const std::string &name,
                 const std::string &string_quadrature,
                 const std::string &string_order)
  :
  ParameterAcceptor(name),
  string_quadrature(string_quadrature),
  string_order(string_order)
{};

template <int dim>
void
ParsedQuadrature<dim>::
declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &string_quadrature,
                "Quadrature to generate", string_quadrature,
                Patterns::Selection("gauss|midpoint|milne|simpson|trapez|weddle"),
                "Description");

  add_parameter(prm, &order,
                "Quadrature order", string_order,
                Patterns::Integer(),
                "Description");
};


D2K_NAMESPACE_CLOSE
