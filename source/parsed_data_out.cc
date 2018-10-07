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

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_tools.h>

#include <deal2lkit/parsed_data_out.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/data_out.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>


D2K_NAMESPACE_OPEN

template <int dim, int spacedim>
ParsedDataOut<dim, spacedim>::ParsedDataOut(
  const std::string & name,
  const std::string & output_format,
  const unsigned int &subdivisions,
  const std::string & incremental_run_prefix,
  const std::string & base_name_input,
  const std::string & files_to_save,
  const MPI_Comm &    comm)
  : ParameterAcceptor(name)
  , comm(comm)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(comm))
  , this_mpi_process(Utilities::MPI::this_mpi_process(comm))
  , output_format(output_format)
  , subdivisions(subdivisions)
  , base_name(base_name_input)
  , incremental_run_prefix(incremental_run_prefix)
  , files_to_save(files_to_save)
{
  initialized = false;
}

template <int dim, int spacedim>
void
ParsedDataOut<dim, spacedim>::declare_parameters(ParameterHandler &prm)
{
  add_parameter(
    prm, &base_name, "Problem base name", base_name, Patterns::Anything());

  add_parameter(prm,
                &incremental_run_prefix,
                "Incremental run prefix",
                incremental_run_prefix,
                Patterns::Anything());

  add_parameter(prm,
                &files_to_save,
                "Files to save in run directory",
                files_to_save,
                Patterns::Anything());

  add_parameter(prm,
                &output_partitioning,
                "Output partitioning",
                "false",
                Patterns::Bool());
  add_parameter(prm,
                &solution_names,
                "Solution names",
                "u",
                Patterns::Anything(),
                "Comma separated list of names for the components. If a "
                "name is repeated, then the repeated names are grouped into "
                "vectors.");

  add_parameter(prm,
                &output_format,
                "Output format",
                output_format,
                Patterns::Selection(DataOutBase::get_output_format_names()));

  add_parameter(prm,
                &subdivisions,
                "Subdivisions",
                std::to_string(subdivisions),
                Patterns::Integer(0));
}

template <int dim, int spacedim>
void
ParsedDataOut<dim, spacedim>::parse_parameters_call_back()
{
  if (incremental_run_prefix != "")
    {
      path_solution_dir =
        get_next_available_directory_name(incremental_run_prefix);
      // The use of the barrier is
      //  to avoid the case of a processor below the master node.
#ifdef DEAL_II_WITH_MPI
      MPI_Barrier(comm);
#endif
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        create_directory(path_solution_dir);
#ifdef DEAL_II_WITH_MPI
      MPI_Barrier(comm);
#endif
      path_solution_dir += "/";
    }
  else
    {
      path_solution_dir = "./";
    }
}

template <int dim, int spacedim>
void
ParsedDataOut<dim, spacedim>::prepare_data_output(
  const DoFHandler<dim, spacedim> &dh,
  const std::string &              suffix)
{
  deallog.push("PrepareOutput");
  data_out = SP(new DataOut<dim, DoFHandler<dim, spacedim>>());
  data_out->set_default_format(DataOutBase::parse_output_format(output_format));

  current_name = path_solution_dir + base_name + suffix;

  std::string fname = current_name;

  if (data_out->default_suffix() != "")
    {
      // If the output is needed and we have many processes, just output
      // the one we need *in intermediate format*.
      if (n_mpi_processes > 1)
        {
          fname += ("." + Utilities::int_to_string(this_mpi_process, 2) + "." +
                    Utilities::int_to_string(n_mpi_processes, 2) +
                    data_out->default_suffix());
        }
      else
        {
          fname += data_out->default_suffix();
        }

      deallog << "Will write on file: " << fname.c_str() << std::endl;
      output_file.open(fname.c_str());
      AssertThrow(output_file, ExcIO());
      data_out->attach_dof_handler(dh);

      if (n_mpi_processes > 1)
        {
          // Output the partitioning
          if (output_partitioning)
            {
              deallog << "Writing partitioning" << std::endl;
              Vector<float> partitioning(
                dh.get_triangulation().n_active_cells());
              for (unsigned int i = 0; i < partitioning.size(); ++i)
                partitioning(i) = this_mpi_process;
              static Vector<float> static_partitioning;
              static_partitioning.swap(partitioning);
              data_out->add_data_vector(static_partitioning, "partitioning");
            }
        }
    }
  initialized = true;
  deallog.pop();
}



template <int dim, int spacedim>
void
ParsedDataOut<dim, spacedim>::write_data_and_clear(
  const Mapping<dim, spacedim> &mapping)
{
  if (output_format == "none")
    return;

  if (files_to_save != "")
    {
      std::vector<std::string> files;
      files = Utilities::split_string_list(files_to_save, '%');
      for (unsigned int i = 0; i < files.size(); ++i)
        {
          if (this_mpi_process == 0)
            copy_files(files[i], path_solution_dir);
        }
    }
  AssertThrow(initialized, ExcNotInitialized());
  AssertThrow(output_file, ExcIO());
  deallog.push("WritingData");
  if (data_out->default_suffix() != "")
    {
      data_out->build_patches(
        mapping,
        subdivisions,
        DataOut<dim, DoFHandler<dim, spacedim>>::curved_inner_cells);
      data_out->write(output_file);
      deallog << "Wrote output file." << std::endl;

      if (this_mpi_process == 0 && n_mpi_processes > 1 &&
          data_out->default_suffix() == ".vtu")
        {
          std::vector<std::string> filenames;
          for (unsigned int i = 0; i < n_mpi_processes; ++i)
            filenames.push_back(current_name + "." +
                                Utilities::int_to_string(i, 2) + "." +
                                Utilities::int_to_string(n_mpi_processes, 2) +
                                data_out->default_suffix());

          std::ofstream master_output((current_name + ".pvtu").c_str());
          data_out->write_pvtu_record(master_output, filenames);
        }
    }
  data_out = 0;
  deallog << "Reset output." << std::endl;
  initialized = false;
  output_file.close();
  deallog.pop();
}


D2K_NAMESPACE_CLOSE

template class deal2lkit::ParsedDataOut<1, 1>;
template class deal2lkit::ParsedDataOut<1, 2>;
template class deal2lkit::ParsedDataOut<1, 3>;
template class deal2lkit::ParsedDataOut<2, 2>;
template class deal2lkit::ParsedDataOut<2, 3>;
template class deal2lkit::ParsedDataOut<3, 3>;
