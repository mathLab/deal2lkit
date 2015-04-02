#ifndef lpcm_post_processor_h
#define lpcm_post_processor_h

#include <fstream>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_vector.h>

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

// #include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/config.h>

#include <deal.II/base/parameter_handler.h>

#include "vector_space.h"
#include "output_processor.h"
#include "error_handler.h"
#include <map>

template <int dim, typename VECTOR=BlockVector<double> >
class PostProcessor : public Subscriptor
{
public:
  /** The constructor takes the mpi initialization stuff. */
  PostProcessor (const unsigned int n_mpi_processes = 1,
                 const unsigned int this_mpi_process = 0);

  /** Initialize the given values for the paramter file. */
  static void declare_parameters(ParameterHandler &prm,
                                 unsigned int ntables=1);

  /** Parse the given parameter handler. */
  void parse_parameters(ParameterHandler &prm);

  /** Calculate the error of the numeric solution in variuous norms. Store
      the result in the given table. */
  void error_from_exact(const VectorSpace<dim> &vspace,
                        const VECTOR &solution,
                        const Function<dim> &exact,
                        unsigned int table_no = 0,
                        double dt=0.);

  /** Difference between two solutions in two different vector spaces. */
  void difference(const VectorSpace<dim> &, const VECTOR &,
                  const VectorSpace<dim> &, const VECTOR &,
                  unsigned int table_no = 0, double dt=0.);

  /** Difference between two solutions in the same vector space. */
  void difference(const VectorSpace<dim> &, const VECTOR &,
                  const VECTOR &, unsigned int table_no = 0, double dt=0.);

  /** Output solution. Using some standard format. This is now deprecated. */
  void output_solution(const VectorSpace<dim> &, const VECTOR &, const std::string &);

  /** Prepare to output data on the given file. This will initialize
      the data_out object and a file with a filename that is the
      combination of the @p base_name, @p filename, eventually a
      processor number and the output suffix.  */
  void prepare_data_output(const DoFHandler<dim> &dh, const std::string &filename);

  /** Add the given vector to the output file. Prior to calling this
      method, you have to call the prepare_data_output method. The
      string can be a comma separated list of components, or a single
      description. In this latter case, a progressive number per
      component is added in the end. */
  void add_data_vector(const VECTOR &data_vector, const std::string &desc);

  /** Actually write the file. Once the data_out has been prepared,
      vectors have been added, the data can be written to a file. This
      is done in this class. At the end of this function call,
      data_out and output_file are in a pristine situation, and the
      process can be started again.*/
  void write_data_and_clear();

  /** Dump vector. Just dump the vector to a file. This is useful to
      read back in the future and restart from where we left off.*/
  void dump_vector(const VECTOR &, const std::string &);

  /** By default output first table. */
  void output_table(const unsigned int table_no=0);

  /** Base name for output files. */
  std::string base_name;

private:
  /** Error results.*/
  std::vector<ConvergenceTable>  tables;

  /** Headers for tables and output. Contains the name of the solution
      components. */
  std::vector<std::string> headers;

  /** Headers for latex tables. Contains the name of the solution
      components. */
  std::vector<std::string> latex_headers;

  /** Captions for latex. */
  std::vector<std::string> latex_captions;

  /** Names of the tables. */
  std::vector<std::string> names;

  /** Type of error to compute per components, separated by ";"
      - 0,0 : don't compute this error;
      - 0,1 : L2 norm.
      - 1,0 : H1 seminorm.
      - 1,1 : H1 norm. Will print also the L2 Norm.
  */
  std::vector<std::vector<NormFlags> > types;

  /** The parameters have been read. */
  bool initialized;

  /** Write the solution. */
  bool write_solution;

  /** Solution format. */
  std::string solution_format;

  /** Compute the error. If this is false, all functions regarding
      errors are disabled and don't do anything.*/
  bool compute_error;

  /** Add convergence rates. */
  std::vector<bool> add_rates;

  /** Write the error files. */
  bool write_error;

  /** Output the error file also on screen. */
  bool output_error;

  /** Output the partitioning of the domain. */
  bool output_partitioning;

  /** The error file format. */
  std::string error_file_format;

  std::vector<std::map<std::string, bool> > extras;

  /** Number of MPI processes. */
  const unsigned int n_mpi_processes;

  /** Id of MPI process. */
  const unsigned int this_mpi_process;

  /** Output file. */
  std::ofstream output_file;

  /** Outputs only the data that refers to this process. */
  FilteredDataOut<dim, DoFHandler<dim> > data_out;
};

#endif
