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

#include <deal2lkit/parsed_function.h>
#include <deal2lkit/utilities.h>

using namespace dealii;

D2K_NAMESPACE_OPEN

template <int dim>
ParsedFunction<dim>::ParsedFunction(const std::string & name,
                                    const unsigned int &n_components,
                                    const std::string & default_exp,
                                    const std::string & default_const)
  : ParameterAcceptor(name)
  , dealii::Functions::ParsedFunction<dim>(n_components)
  , default_exp(default_exp)
  , default_const(default_const)
  , n_components(n_components)
{}


template <int dim>
void
ParsedFunction<dim>::declare_parameters(dealii::ParameterHandler &prm)
{
  dealii::Functions::ParsedFunction<dim>::declare_parameters(prm, n_components);
  if (default_exp != "")
    prm.set("Function expression", default_exp);
  if (default_const != "")
    prm.set("Function constants", default_const);
}


template <int dim>
void
ParsedFunction<dim>::parse_parameters(dealii::ParameterHandler &prm)
{
  dealii::Functions::ParsedFunction<dim>::parse_parameters(prm);
}


template class ParsedFunction<1>;
template class ParsedFunction<2>;
template class ParsedFunction<3>;

D2K_NAMESPACE_CLOSE
