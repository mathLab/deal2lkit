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
ParsedQuadrature(const std::string    &name,
                 const std::string    &quadrature_type,
                 const unsigned int   order)
  :
  ParameterAcceptor(name),
  quadrature_type(quadrature_type),
  order(order)
{};

template <int dim>
void
ParsedQuadrature<dim>::
declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &quadrature_type,
                "Quadrature to generate", quadrature_type,
                Patterns::Selection("gauss|midpoint|milne|simpson|trapez|weddle"),
                "Description");

  add_parameter(prm, &order,
                "Quadrature order", std::to_string(order),
                Patterns::Integer(0),
                "Description");
};

template <int dim>
void
ParsedQuadrature<dim>::
parse_parameters(ParameterHandler &prm)
{};

template <int dim>
void
ParsedQuadrature<dim>::
parse_parameters_call_back()
{};

template <int dim>
Quadrature<dim>
ParsedQuadrature<dim>::
get_quadrature()
{
  return QuadratureSelector<dim>(quadrature_type, order);
};


D2K_NAMESPACE_CLOSE

template class deal2lkit::ParsedQuadrature<1>;
template class deal2lkit::ParsedQuadrature<2>;
template class deal2lkit::ParsedQuadrature<3>;
