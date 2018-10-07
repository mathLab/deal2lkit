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

#include <deal.II/base/quadrature.h>

#include <deal2lkit/parameter_acceptor.h>

using namespace dealii;


D2K_NAMESPACE_OPEN

/**
 * A deal2lkit wrapper for an iterated version of
 * dealii::QuadratureSelector.  The template integer specifies the
 * dimension of the quadrature.
 */
template <int dim>
class ParsedQuadrature : public Quadrature<dim>, public ParameterAcceptor
{
public:
  /**
   * Constructor: takes an optional name for the section.
   * Moreover, it Takes the name of the quadrature rule @p quadrature_type that
   * could be chosen among following names:
   * -  gauss;
   * -  midpoint;
   * -  milne;
   * -  simpson;
   * -  trapez;
   * -  weddle;
   * If @p quadrature_type is "gauss" it is possible to specify the number of
   * quadrature points in each coordinate direction (@p order).
   *
   * The last parameter specify how many times this quadrature should
   * be repeated, i.e., in one space dimension, the given base formula
   * is copied and scaled onto a given number of subintervals of
   * length 1/@p repetitions. In more than one space dimension, the
   * resulting quadrature formula is constructed by tensor product of
   * the respective iterated quadrature formula in one space
   * dimension.
   */
  ParsedQuadrature(const std::string &name            = "",
                   const std::string &quadrature_type = "gauss",
                   const unsigned int order           = 3,
                   const unsigned int repetitions     = 1);

  /**
   * Declare quadrature type and quadrature options.
   */
  virtual void declare_parameters(ParameterHandler &prm);

  /**
   * Fill information about blocks after parsing the parameters.
   */
  virtual void parse_parameters_call_back();

private:
  /**
   * Name of the quadrature of the quadrature rule: "gauss", "midpoint",
   * "milne", "simpson", "trapez", or "weddle".
   */
  std::string quadrature_type;

  /**
   * In one space dimension, the given base formula is copied and scaled onto a
   * given number of subintervals of length 1/@p repetitions. In more than one
   * space dimension, the resulting quadrature formula is constructed in the
   * usual way by building the tensor product of the respective iterated
   * quadrature formula in one space dimension.
   */
  unsigned int repetitions;

  /**
   * Number of quadrature points in each coordinate direction.
   * This variable is only valid for gauss rule.
   */
  unsigned int order;
};

D2K_NAMESPACE_CLOSE


#endif
