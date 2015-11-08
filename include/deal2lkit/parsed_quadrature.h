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

#ifndef _d2k_parsed_quadrature
#define _d2k_parsed_quadrature

#include <deal.II/base/quadrature_selector.h>

#include <deal2lkit/parameter_acceptor.h>

using namespace dealii;


D2K_NAMESPACE_OPEN

/**
 * TODO:
 */
template<int dim>
class ParsedQuadrature : public ParameterAcceptor
{
public:
  /**
   * TODO:
   */
  ParsedQuadrature( const std::string &name="",
                    const std::string &string_quadrature="gauss",
                    const std::string &string_order="3");

  /**
   * Declare quadrature type and solver options.
   */
  virtual void declare_parameters(ParameterHandler &prm);

  /**
   * Parse quadrature type and solver options.
   */
  virtual void parse_parameters(ParameterHandler &prm);

  /**
   * Initialize internal variables.
   */
  virtual void parse_parameters_call_back();

private:
  std::string string_quadrature;
  std::string string_order;

  unsigned int order;
};

D2K_NAMESPACE_CLOSE


#endif
