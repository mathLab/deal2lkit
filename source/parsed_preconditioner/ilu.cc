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

using namespace dealii;

D2K_NAMESPACE_OPEN

ParsedILUPreconditioner::ParsedILUPreconditioner(const std::string & name,
                                                 const unsigned int &ilu_fill,
                                                 const double &      ilu_atol,
                                                 const double &      ilu_rtol,
                                                 const unsigned int &overlap)
  : ParameterAcceptor(name)
  , PreconditionILU()
  , ilu_fill(ilu_fill)
  , ilu_atol(ilu_atol)
  , ilu_rtol(ilu_rtol)
  , overlap(overlap)
{}

void
ParsedILUPreconditioner::declare_parameters(ParameterHandler &prm)
{
  prm.add_parameter("Fill-in", ilu_fill, "Additional fill-in.");

  prm.add_parameter("ILU atol",
                    ilu_atol,
                    "The amount of perturbation to add to diagonal entries.");

  prm.add_parameter("ILU rtol",
                    ilu_rtol,
                    "Scaling factor for diagonal entries.");

  prm.add_parameter("Overlap", overlap, "Overlap between processors.");
}

template <typename Matrix>
void
ParsedILUPreconditioner::initialize_preconditioner(const Matrix &matrix)
{
  TrilinosWrappers::PreconditionILU::AdditionalData data;

  data.ilu_fill = ilu_fill;
  data.ilu_atol = ilu_atol;
  data.ilu_rtol = ilu_rtol;
  data.overlap  = overlap;
  this->initialize(matrix, data);
}
D2K_NAMESPACE_CLOSE

template void
deal2lkit::ParsedILUPreconditioner::initialize_preconditioner<
  dealii::TrilinosWrappers::SparseMatrix>(
  const dealii::TrilinosWrappers::SparseMatrix &);

#endif
