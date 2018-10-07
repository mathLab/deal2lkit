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

#include <deal2lkit/parsed_preconditioner/jacobi.h>

#ifdef DEAL_II_WITH_TRILINOS

D2K_NAMESPACE_OPEN

ParsedJacobiPreconditioner::ParsedJacobiPreconditioner(
  const std::string & name,
  const double &      omega,
  const double &      min_diagonal,
  const unsigned int &n_sweeps) :
  ParameterAcceptor(name),
  PreconditionJacobi(),
  omega(omega),
  min_diagonal(min_diagonal),
  n_sweeps(n_sweeps)
{}

void ParsedJacobiPreconditioner::declare_parameters(ParameterHandler &prm)
{
  add_parameter(
    prm,
    &omega,
    "Omega",
    std::to_string(omega),
    Patterns::Double(0.0),
    "This specifies the relaxation parameter in the Jacobi preconditioner.");
  add_parameter(
    prm,
    &min_diagonal,
    "Min Diagonal",
    std::to_string(min_diagonal),
    Patterns::Double(0.0),
    "This specifies the minimum value the diagonal elements should\n"
    "have. This might be necessary when the Jacobi preconditioner is used\n"
    "on matrices with zero diagonal elements. In that case, a straight-\n"
    "forward application would not be possible since we would divide by\n"
    "zero.");
  add_parameter(
    prm,
    &n_sweeps,
    "Number of sweeps",
    std::to_string(n_sweeps),
    Patterns::Integer(0),
    "Sets how many times the given operation should be applied during the\n"
    "vmult() operation.");
}

template <typename Matrix>
void ParsedJacobiPreconditioner::initialize_preconditioner(const Matrix &matrix)
{
  TrilinosWrappers::PreconditionJacobi::AdditionalData data;

  data.omega        = omega;
  data.min_diagonal = min_diagonal;
  data.n_sweeps     = n_sweeps;
  this->initialize(matrix, data);
}
D2K_NAMESPACE_CLOSE

template void deal2lkit::ParsedJacobiPreconditioner::initialize_preconditioner<
  dealii::TrilinosWrappers::SparseMatrix>(
  const dealii::TrilinosWrappers::SparseMatrix &);

#endif
