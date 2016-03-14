
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

#include <deal2lkit/parsed_preconditioner_amg.h>

#ifdef DEAL_II_WITH_TRILINOS


#include <deal.II/dofs/dof_tools.h>


D2K_NAMESPACE_OPEN

ParsedAMGPreconditioner::ParsedAMGPreconditioner( const std::string &name,
                                                  const bool &elliptic,
                                                  const bool &higher_order_elements,
                                                  const unsigned int &n_cycles,
                                                  const bool &w_cycle,
                                                  const double &aggregation_threshold,
                                                  const std::string &var_const_modes,
                                                  const unsigned int &smoother_sweeps,
                                                  const unsigned int &smoother_overlap,
                                                  const bool &output_details,
                                                  const std::string &smoother_type,
                                                  const std::string &coarse_type
                                                ) :
  ParameterAcceptor(name),
  PreconditionAMG(),
  elliptic(elliptic),
  higher_order_elements(higher_order_elements),
  n_cycles(n_cycles),
  w_cycle(w_cycle),
  aggregation_threshold(aggregation_threshold),
  var_const_modes(var_const_modes),
  smoother_sweeps(smoother_sweeps),
  smoother_overlap(smoother_overlap),
  output_details(output_details),
  smoother_type(smoother_type),
  coarse_type(coarse_type)
{}

void ParsedAMGPreconditioner::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &elliptic, "Elliptic", elliptic ? "true":"false",
                Patterns::Bool(),
                "Determines whether the AMG preconditioner should be optimized for\n"
                "elliptic problems (ML option smoothed aggregation SA, using a\n"
                "Chebyshev smoother) or for non-elliptic problems (ML option\n"
                "non-symmetric smoothed aggregation NSSA, smoother is SSOR with\n"
                "underrelaxation");

  add_parameter(prm, &higher_order_elements, "High Order Elements", higher_order_elements ? "true":"false",
                Patterns::Bool(),
                "Determines whether the matrix that the preconditioner is built upon is\n"
                "generated from linear or higher-order elements.");

  add_parameter(prm, &n_cycles, "Number of cycles", std::to_string(n_cycles),
                Patterns::Integer(0),
                "Defines how many multigrid cycles should be performed by the\n"
                "preconditioner.");

  add_parameter(prm, &w_cycle, "w-cycle", w_cycle ? "true":"false",
                Patterns::Bool(),
                "efines whether a w-cycle should be used instead of the standard\n"
                "setting of a v-cycle.");

  add_parameter(prm, &aggregation_threshold, "Aggregation threshold", std::to_string(aggregation_threshold),
                Patterns::Double(0.0),
                "This threshold tells the AMG setup how the coarsening should be\n"
                "performed. In the AMG used by ML, all points that strongly couple with\n"
                "the tentative coarse-level point form one aggregate. The term strong\n"
                "coupling is controlled by the variable aggregation_threshold, meaning\n"
                "that all elements that are not smaller than aggregation_threshold\n"
                "times the diagonal element do couple strongly.");

  add_parameter(prm, &var_const_modes,
                "Variable related to constant modes", var_const_modes,
                Patterns::Anything(),
                "Specifies the variable associated to the constant modes (near null\n"
                "space) of the matrix. In the case @p var_const_modes\n"
                "is equal to \"none\", constant modes will not be\n"
                "computed.");

  add_parameter(prm, &smoother_sweeps, "Smoother sweeps", std::to_string(smoother_sweeps),
                Patterns::Integer(0),
                "Determines how many sweeps of the smoother should be performed. When\n"
                "the flag elliptic is set to true, i.e., for elliptic or almost\n"
                "elliptic problems, the polynomial degree of the Chebyshev smoother is\n"
                "set to smoother_sweeps. The term sweeps refers to the number of\n"
                "matrix-vector products performed in the Chebyshev case. In the\n"
                "non-elliptic case, smoother_sweeps sets the number of SSOR relaxation\n"
                "sweeps for post-smoothing to be performed.");

  add_parameter(prm, &smoother_overlap, "Smoother overlap", std::to_string(smoother_overlap),
                Patterns::Integer(0),
                "Determines the overlap in the SSOR/Chebyshev error smoother when run\n"
                "in parallel.");

  add_parameter(prm, &output_details, "Output details", output_details ? "true":"false",
                Patterns::Bool(),
                "If this flag is set to true, then internal information from the ML\n"
                "preconditioner is printed to screen. This can be useful when debugging\n"
                "the preconditioner.");

  add_parameter(prm, &smoother_type, "Smoother type", smoother_type,
                Patterns::Selection("|Aztec|IFPACK|Jacobi"
                                    "|ML symmetric Gauss-Seidel|symmetric Gauss-Seidel"
                                    "|ML Gauss-Seidel|Gauss-Seidel|block Gauss-Seidel"
                                    "|symmetric block Gauss-Seidel|Chebyshev|MLS|Hiptmair"
                                    "|Amesos-KLU|Amesos-Superlu|Amesos-UMFPACK|Amesos-Superludist"
                                    "|Amesos-MUMPS|user-defined|SuperLU|IFPACK-Chebyshev|self"
                                    "|do-nothing|IC|ICT|ILU|ILUT|Block Chebyshev"
                                    "|IFPACK-Block Chebyshev"),
                "Determines which smoother to use for the AMG cycle.");

  add_parameter(prm, &coarse_type, "Coarse type", coarse_type,
                Patterns::Selection("|Aztec|IFPACK|Jacobi"
                                    "|ML symmetric Gauss-Seidel|symmetric Gauss-Seidel"
                                    "|ML Gauss-Seidel|Gauss-Seidel|block Gauss-Seidel"
                                    "|symmetric block Gauss-Seidel|Chebyshev|MLS|Hiptmair"
                                    "|Amesos-KLU|Amesos-Superlu|Amesos-UMFPACK|Amesos-Superludist"
                                    "|Amesos-MUMPS|user-defined|SuperLU|IFPACK-Chebyshev|self"
                                    "|do-nothing|IC|ICT|ILU|ILUT|Block Chebyshev"
                                    "|IFPACK-Block Chebyshev"),
                "Determines which solver to use on the coarsest level. The same\n"
                "settings as for the smoother type are possible.");
}

template<typename Matrix>
void ParsedAMGPreconditioner::initialize_preconditioner( const Matrix &matrix)
{
  TrilinosWrappers::PreconditionAMG::AdditionalData data;

  data.elliptic = elliptic;
  data.higher_order_elements = higher_order_elements;
  data.n_cycles = n_cycles;
  data.w_cycle = w_cycle;
  data.aggregation_threshold = aggregation_threshold;
  data.constant_modes = std::vector<std::vector<bool> > (0);
  data.smoother_sweeps = smoother_sweeps;
  data.smoother_overlap = smoother_overlap;
  data.output_details = output_details;
  data.smoother_type = smoother_type.c_str();
  data.coarse_type = coarse_type.c_str();
  this->initialize(matrix, data);
}

template<int dim, int spacedim, typename Matrix>
void ParsedAMGPreconditioner::initialize_preconditioner( const Matrix &matrix,
                                                         const ParsedFiniteElement<dim, spacedim> &fe,
                                                         const DoFHandler<dim, spacedim> &dh)
{
  TrilinosWrappers::PreconditionAMG::AdditionalData data;

  data.elliptic = elliptic;
  data.higher_order_elements = higher_order_elements;
  data.n_cycles = n_cycles;
  data.w_cycle = w_cycle;
  data.aggregation_threshold = aggregation_threshold;
  if (var_const_modes == "none")
    {
      data.constant_modes = std::vector<std::vector<bool> > (0);
    }
  else
    {
      std::vector< std::vector< bool > >  constant_modes;
      unsigned int pos = fe.get_first_occurence(var_const_modes);
      bool is_vec = fe.is_vector(var_const_modes);
      if (is_vec)
        {
          FEValuesExtractors::Vector components(pos);
          DoFTools::extract_constant_modes<DoFHandler<dim, spacedim> > (dh, dh.get_fe().component_mask(components),
              constant_modes);
        }
      else
        {
          FEValuesExtractors::Scalar components(pos);
          DoFTools::extract_constant_modes<DoFHandler<dim, spacedim> > (dh, dh.get_fe().component_mask(components),
              constant_modes);
        }
      data.constant_modes = constant_modes;
    }
  data.smoother_sweeps = smoother_sweeps;
  data.smoother_overlap = smoother_overlap;
  data.output_details = output_details;
  data.smoother_type = smoother_type.c_str();
  data.coarse_type = coarse_type.c_str();
  this->initialize(matrix, data);
}

D2K_NAMESPACE_CLOSE

template void deal2lkit::ParsedAMGPreconditioner::initialize_preconditioner<1,1,dealii::TrilinosWrappers::SparseMatrix>(
  const dealii::TrilinosWrappers::SparseMatrix &,
  const ParsedFiniteElement<1,1> &,
  const DoFHandler<1,1> &);
template void deal2lkit::ParsedAMGPreconditioner::initialize_preconditioner<2,2,dealii::TrilinosWrappers::SparseMatrix>(
  const dealii::TrilinosWrappers::SparseMatrix &,
  const ParsedFiniteElement<2,2> &,
  const DoFHandler<2,2> &);
template void deal2lkit::ParsedAMGPreconditioner::initialize_preconditioner<3,3,dealii::TrilinosWrappers::SparseMatrix>(
  const dealii::TrilinosWrappers::SparseMatrix &,
  const ParsedFiniteElement<3,3> &,
  const DoFHandler<3,3> &);

template void deal2lkit::ParsedAMGPreconditioner::initialize_preconditioner<dealii::TrilinosWrappers::SparseMatrix>(
  const dealii::TrilinosWrappers::SparseMatrix &);

#endif
