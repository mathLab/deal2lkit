#ifndef __dealii_sak_parsed_data_out_h
#define __dealii_sak_parsed_data_out_h

#include <fstream>

#include <deal.II/base/config.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>

#include "parameter_acceptor.h"


template <int dim, int spacedim=dim>
class ParsedDataOut : public ParameterAcceptor
{
public:
  /** Optional name for parameter section. */
  ParsedDataOut (const std::string &name="",
                 const std::string &default_format="vtu",
                 const MPI_Comm &comm=MPI_COMM_WORLD);

  /** Initialize the given values for the paramter file. */
  virtual void declare_parameters(ParameterHandler &prm);

  /** Parse the given parameters. */
  virtual void parse_parameters(ParameterHandler &prm);

  /** Prepare to output data on the given file. This will initialize
      the data_out object and a file with a filename that is the
      combination of the @p base_name, the optional @p suffix,
      eventually a processor number and the output suffix.  */
  void prepare_data_output(const DoFHandler<dim, spacedim> &dh,
                           const std::string &suffix="");

  /** Add the given vector to the output file. Prior to calling this
      method, you have to call the prepare_data_output method. The
      string can be a comma separated list of components, or a single
      description. In this latter case, a progressive number per
      component is added in the end. */
  template<typename VECTOR>
  void add_data_vector(const VECTOR &data_vector, const std::string &desc);

  /** Actually write the file. Once the data_out has been prepared,
      vectors have been added, the data can be written to a file. This
      is done in this class. At the end of this function call,
      data_out and output_file are in a pristine situation, and the
      process can be started again.*/
  void write_data_and_clear(const Mapping<dim,spacedim> &mapping=StaticMappingQ1<dim,spacedim>::mapping);

private:
  /** Initialization flag.*/
  bool initialized;

  /** MPI communicator. */
  const MPI_Comm &comm;

  /** Number of processes. */
  const unsigned int n_mpi_processes;

  /** My mpi process. */
  const unsigned int this_mpi_process;

  /** Default format at construction time. */
  const std::string default_format;

  /** Base name for output files. This base is used to generate all
      filenames. */
  std::string base_name;

  /** Solution names. */
  std::string solution_names;

  /** Current output name. When preparing data_out, this name will
      contain the base for the current output. This allows the user to
      use a different output name in different part of the program. */
  std::string current_name;

  /** Output the partitioning of the domain. */
  bool output_partitioning;

  /** Output file. */
  std::ofstream output_file;

  /** Outputs only the data that refers to this process. */
  DataOut<dim, DoFHandler<dim, spacedim> > data_out;
};


// ============================================================
// Template specializations
// ============================================================

template <int dim, int spacedim>
template<typename VECTOR>
void ParsedDataOut<dim,spacedim>::add_data_vector(const VECTOR &data_vector,
                                                  const std::string &desc)
{
  AssertThrow(initialized, ExcNotInitialized());
  deallog.push("AddingData");
  std::vector<std::string> dd = Utilities::split_string_list(desc);
  if (data_out.default_suffix() != "")
    {
      if (dd.size() ==1 )
        data_out.add_data_vector (data_vector, desc);
      else
        data_out.add_data_vector (data_vector, dd);
      deallog << "Added data: " << desc << std::endl;
    }
  deallog.pop();
}



#endif
