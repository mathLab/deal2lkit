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

#ifndef d2k_parsed_function_h
#define d2k_parsed_function_h

#include <deal.II/base/parsed_function.h>

#include <deal2lkit/config.h>
#include <deal2lkit/parameter_acceptor.h>

using namespace dealii;

D2K_NAMESPACE_OPEN

/**
 * A deal2lkit wrapper for dealii::Functions::ParsedFunction. The
 * template integer specify the dimension of points this function
 * accepts.
 */
template <int dim>
class ParsedFunction : public ParameterAcceptor,
                       public Functions::ParsedFunction<dim>
{
public:
  /**
   * Constructor: takes an optional name for the section. If the
   * optional expression string is given, than it is used to set the
   * expression as soon as the parameters are declared.
   *
   * The n_components parameter is used to determine the number of
   * components of this function. It defaults to one, and once it is
   * set in the constructor, it cannot be changed by parameter file.
   * The default_expression, instead, can be changed by acting on the
   * parameter file. An exception will be thrown if the number of
   * components of the expression does not coincide with the parameter
   * passed at construction time to this function.
   */
  ParsedFunction(const std::string & name          = "",
                 const unsigned int &n_components  = 1,
                 const std::string & default_exp   = "",
                 const std::string & default_const = "");

  /**
   * Calls the underlying function of ParsedFunction.
   */
  virtual void
  declare_parameters(ParameterHandler &prm);

  /**
   * Calls the underlying function of ParsedFunction.
   */
  virtual void
  parse_parameters(ParameterHandler &prm);


private:
  /**
   * Default expression of this function. "
   */
  const std::string  default_exp;
  const std::string  default_const;
  const unsigned int n_components;
};

D2K_NAMESPACE_CLOSE

#endif
