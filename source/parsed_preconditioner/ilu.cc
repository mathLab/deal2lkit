//-----------------------------------------------------------
//
//    Copyright (C) 2016 by the deal2lkit authors
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

#include <deal2lkit/parsed_preconditioner/ilu.h>

#ifdef DEAL_II_WITH_TRILINOS

D2K_NAMESPACE_OPEN

ParsedPreconditionerILU::ParsedPreconditionerILU(const std::string &name,
                                                 const unsigned int &ilu_fill,
                                                 const double &ilu_atol,
                                                 const double &ilu_rtol,
                                                 const unsigned int &overlap
                                                ):
  ParameterAcceptor(name),
  PreconditionILU(),
  ilu_fill(ilu_fill),
  ilu_atol(ilu_atol),
  ilu_rtol(ilu_rtol),
  overlap(overlap)
{}

void ParsedPreconditionerILU::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &ilu_fill, "Fill-in", std::to_string(ilu_fill),
                Patterns::Integer(0),
                "Additional fill-in.");
  add_parameter(prm, &ilu_atol, "ILU atol", std::to_string(ilu_atol),
                Patterns::Double(0.0),
                "The amount of perturbation to add to diagonal entries.");
  add_parameter(prm, &ilu_rtol, "ILU rtol", std::to_string(ilu_rtol),
                Patterns::Double(0.0),
                "Scaling factor for diagonal entries.");
  add_parameter(prm, &overlap, "Overlap", std::to_string(overlap),
                Patterns::Integer(0),
                "Overlap between processors.");
}

template<typename Matrix>
void ParsedPreconditionerILU::initialize_preconditioner( const Matrix &matrix)
{
  TrilinosWrappers::PreconditionILU::AdditionalData data;

  data.ilu_fill = ilu_fill;
  data.ilu_atol = ilu_atol;
  data.ilu_rtol = ilu_rtol;
  data.overlap = overlap;
  this->initialize(matrix, data);
}
D2K_NAMESPACE_CLOSE

template void deal2lkit::ParsedPreconditionerILU::initialize_preconditioner<dealii::TrilinosWrappers::SparseMatrix>(
  const dealii::TrilinosWrappers::SparseMatrix &);

#endif
