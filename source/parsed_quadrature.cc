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
#include <deal.II/base/quadrature_selector.h>

D2K_NAMESPACE_OPEN

template <int dim>
ParsedQuadrature<dim>::
ParsedQuadrature(const std::string    &name,
                 const std::string    &quadrature_type,
                 const unsigned int   order,
                 const unsigned int   repetitions)
  :
  ParameterAcceptor(name),
  quadrature_type(quadrature_type),
  repetitions(repetitions),
  order(order)
{};

template <int dim>
void
ParsedQuadrature<dim>::
declare_parameters(ParameterHandler &prm)
{
  std::string doc_quadrature_type = QuadratureSelector<dim>::get_quadrature_names();

  add_parameter(prm, &quadrature_type,
                "Quadrature to generate", quadrature_type,
                Patterns::Selection(QuadratureSelector<dim>::get_quadrature_names()),
                "Quadrature rule:"+
                doc_quadrature_type
               );

  add_parameter(prm, &order,
                "Quadrature order", std::to_string(order),
                Patterns::Integer(0),
                "The number of quadrature points in each coordinate direction. (Avaible only for gauss otherwise it should be 0)");

  add_parameter(prm, &repetitions,
                "Number of repetitions", std::to_string(repetitions),
                Patterns::Integer(1),
                "In one space dimension, the given base formula is copied and scaled onto a given number of subintervals of length 1/repetitions. In more than one space dimension, the resulting quadrature formula is constructed in the usual way by building the tensor product of the respective iterated quadrature formula in one space dimension.");
};

template <int dim>
void
ParsedQuadrature<dim>::
parse_parameters_call_back()
{
  (Quadrature<dim> &)(* this) =  QIterated<dim>( QuadratureSelector<1>(quadrature_type,order), repetitions );
}

D2K_NAMESPACE_CLOSE

template class deal2lkit::ParsedQuadrature<1>;
template class deal2lkit::ParsedQuadrature<2>;
template class deal2lkit::ParsedQuadrature<3>;
