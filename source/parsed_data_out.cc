#include "../include/parsed_data_out.h"

#include "../include/utilities.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>


#include <deal.II/base/utilities.h>

template <int dim, int spacedim>
ParsedDataOut<dim,spacedim>::ParsedDataOut (const std::string &name,
                                            const std::string &default_format,
                                            const std::string &incremental_run_prefix,
                                            const std::string &base_name_input,
                                            const MPI_Comm &comm) :
  ParameterAcceptor(name),
  comm(comm),
  n_mpi_processes(Utilities::MPI::n_mpi_processes(comm)),
  this_mpi_process(Utilities::MPI::this_mpi_process(comm)),
  default_format(default_format),
  base_name(base_name_input),
  incremental_run_prefix(incremental_run_prefix)
{
  initialized = false;
}

template <int dim, int spacedim>
void ParsedDataOut<dim,spacedim>::declare_parameters (ParameterHandler &prm)
{
  add_parameter(prm, &base_name, "Problem base name", base_name, Patterns::Anything());
  add_parameter(prm, &incremental_run_prefix, "Problem run dir name", incremental_run_prefix, Patterns::Anything());

  add_parameter(prm, &output_partitioning, "Output partitioning", "false", Patterns::Bool());
  add_parameter(prm, &solution_names, "Solution names", "u", Patterns::Anything(),
                "Comma separated list of names for the components. If a "
                "name is repeated, then the repeated names are grouped into "
                "vectors.");

  prm.enter_subsection("Solution output format");
  DataOut<dim, DoFHandler<dim,spacedim> >::declare_parameters(prm);
  prm.set("Output format", default_format);
  prm.leave_subsection();
}


template <int dim, int spacedim>
void ParsedDataOut<dim,spacedim>::parse_parameters (ParameterHandler &prm)
{
  ParameterAcceptor::parse_parameters(prm);

  prm.enter_subsection("Solution output format");
  data_out.parse_parameters(prm);
  prm.leave_subsection();

  initialized = true;
}

template <int dim, int spacedim>
void ParsedDataOut<dim,spacedim>::parse_parameters_call_back()
{
  if ( incremental_run_prefix != "" )
    {
      path_solution_dir = get_next_available_directory_name(incremental_run_prefix);
      // The use of the barrier is
      //  to avoid the case of a processor below the master node.
#ifdef DEAL_II_WITH_MPI
      MPI_Barrier(comm);
#endif
      if ( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
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
void ParsedDataOut<dim,spacedim>::prepare_data_output(const DoFHandler<dim,spacedim> &dh,
                                                      const std::string &suffix)
{
  AssertThrow(initialized, ExcNotInitialized());
  deallog.push("PrepareOutput");

  current_name = path_solution_dir + base_name + suffix;

  std::string fname = current_name;

  if (data_out.default_suffix() != "")
    {

      // If the output is needed and we have many processes, just output
      // the one we need *in intermediate format*.
      if (n_mpi_processes > 1)
        {
          fname += ("." + Utilities::int_to_string(this_mpi_process, 2) +
                    "." + Utilities::int_to_string(n_mpi_processes, 2) +
                    data_out.default_suffix()) ;
        }
      else
        {
          fname += data_out.default_suffix();
        }

      deallog << "Will write on file: " << fname.c_str() << std::endl;
      output_file.open(fname.c_str());
      AssertThrow(output_file, ExcIO());
      data_out.attach_dof_handler (dh);

      if (n_mpi_processes > 1)
        {
          // Output the partitioning
          if (output_partitioning)
            {
              deallog << "Writing partitioning" << std::endl;
              Vector<float> partitioning(dh.get_tria().n_active_cells());
              for (unsigned int i=0; i<partitioning.size(); ++i)
                partitioning(i) = this_mpi_process;
              static Vector<float> static_partitioning;
              static_partitioning.swap(partitioning);
              data_out.add_data_vector (static_partitioning, "partitioning");
            }
        }
    }

  deallog.pop();
}



template <int dim, int spacedim>
void ParsedDataOut<dim,spacedim>::write_data_and_clear( const std::string &used_files,
                                                        const Mapping<dim,spacedim> &mapping)
{
  copy_files(used_files, path_solution_dir);

  AssertThrow(initialized, ExcNotInitialized());
  AssertThrow(output_file, ExcIO());
  deallog.push("WritingData");
  if (data_out.default_suffix() != "")
    {
      data_out.build_patches(mapping, 0,
                             DataOut<dim, DoFHandler<dim,spacedim> >::curved_inner_cells);
      data_out.write(output_file);
      deallog << "Wrote output file." << std::endl;

      if (this_mpi_process == 0 && n_mpi_processes > 0 && data_out.default_suffix() == ".vtu")
        {
          std::vector<std::string> filenames;
          for (unsigned int i=0; i<n_mpi_processes; ++i)
            filenames.push_back (current_name +
                                 "." + Utilities::int_to_string(i, 2) +
                                 "." + Utilities::int_to_string(n_mpi_processes, 2) +
                                 data_out.default_suffix());

          std::ofstream master_output ((current_name + ".pvtu").c_str());
          data_out.write_pvtu_record (master_output, filenames);
        }

      data_out.clear();
      output_file.close();
      deallog << "Reset output." << std::endl;
    }
  deallog.pop();
}

template class ParsedDataOut<1,1>;
template class ParsedDataOut<1,2>;
template class ParsedDataOut<1,3>;
template class ParsedDataOut<2,2>;
template class ParsedDataOut<2,3>;
template class ParsedDataOut<3,3>;
