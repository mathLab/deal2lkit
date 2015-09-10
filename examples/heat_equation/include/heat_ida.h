#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/base/timer.h>
#include <stdio.h>
#include <stdlib.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/solver_cg.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_math.h>

#include <mpi.h>


#include <nvector/nvector_parallel.h>

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/distributed/solution_transfer.h>

#include <mpi.h>

#include "sundials_interface.h"
#include "dae_time_integrator.h"
#include "parameter_acceptor.h"
#include "parsed_grid_generator.h"
#include "parsed_finite_element.h"
#include "error_handler.h"
#include "parsed_function.h"
#include "parsed_data_out.h"
#include "parsed_dirichlet_bcs.h"

#include <fstream>

typedef typename TrilinosWrappers::MPI::Vector VEC;
template<int dim>
class Heat : public SundialsInterface<VEC>, public ParameterAcceptor
{
public:

  Heat (const MPI_Comm &comm);

  virtual void declare_parameters (ParameterHandler &prm);

  void run ();

  /*********************************************************
   * Public interface from SundialsInterface
   *********************************************************/
  virtual shared_ptr<VEC>
  create_new_vector() const;

  /** Returns the number of degrees of freedom. Pure virtual function. */
  virtual unsigned int n_dofs() const;

  /** This function is called at the end of each iteration step for
   * the ode solver. Once again, the conversion between pointers and
   * other forms of vectors need to be done inside the inheriting
   * class. */
  virtual void output_step(const double t,
                           const VEC &solution,
                           const VEC &solution_dot,
                           const unsigned int step_number,
                           const double h);

  /** This function will check the behaviour of the solution. If it
   * is converged or if it is becoming unstable the time integrator
   * will be stopped. If the convergence is not achived the
   * calculation will be continued. If necessary, it can also reset
   * the time stepper. */
  virtual bool solver_should_restart(const double t,
                                     const VEC &solution,
                                     const VEC &solution_dot,
                                     const unsigned int step_number,
                                     const double h);

  /** For dae problems, we need a
   residual function. */
  virtual int residual(const double t,
                       const VEC &src_yy,
                       const VEC &src_yp,
                       VEC &dst);

  /** Setup Jacobian system and preconditioner. */
  virtual int setup_jacobian(const double t,
                             const VEC &src_yy,
                             const VEC &src_yp,
                             const VEC &residual,
                             const double alpha);


  /** Inverse of the Jacobian vector product. */
  virtual int solve_jacobian_system(const double t,
                                    const VEC &y,
                                    const VEC &y_dot,
                                    const VEC &residual,
                                    const double alpha,
                                    const VEC &src,
                                    VEC &dst) const;



  /** And an identification of the
   differential components. This
   has to be 1 if the
   corresponding variable is a
   differential component, zero
   otherwise.  */
  virtual VEC &differential_components() const;

private:
  void make_grid_fe();
  void setup_dofs (const bool &first_run=true);

  void assemble_jacobian_matrix (const double t,
                                 const VEC &y,
                                 const VEC &y_dot,
                                 const double alpha);

  void refine_mesh ();
  void process_solution ();

  void set_constrained_dofs_to_zero(VEC &v) const;

  const MPI_Comm &comm;

  unsigned int n_cycles;
  unsigned int current_cycle;
  unsigned int initial_global_refinement;
  unsigned int max_time_iterations;
  double fixed_alpha;

  std::string timer_file_name;

  ConditionalOStream        pcout;
  std::ofstream         timer_outfile;
  ConditionalOStream        tcout;

  shared_ptr<Mapping<dim,dim> >             mapping;

  shared_ptr<parallel::distributed::Triangulation<dim,dim> > triangulation;
  shared_ptr<FiniteElement<dim,dim> >       fe;
  shared_ptr<DoFHandler<dim,dim> >          dof_handler;

  ConstraintMatrix                          constraints;

  TrilinosWrappers::SparsityPattern       jacobian_matrix_sp;
  TrilinosWrappers::SparseMatrix          jacobian_matrix;

  TrilinosWrappers::PreconditionAMG       preconditioner;

  VEC        solution;
  VEC        solution_dot;

  mutable VEC        distributed_solution;
  mutable VEC        distributed_solution_dot;


  mutable TimerOutput     computing_timer;


  //ErrorHandler<1>       eh;
  ParsedGridGenerator<dim,dim>   pgg;
  ParsedFiniteElement<dim,dim> fe_builder;

  ParsedFunction<dim, 1>        exact_solution;
  ParsedFunction<dim, 1>        forcing_term;

  ParsedFunction<dim, 1>        initial_solution;
  ParsedFunction<dim, 1>        initial_solution_dot;
  ParsedDirichletBCs<dim,dim,1> dirichlet_bcs;

  ParsedDataOut<dim, dim>                  data_out;

  DAETimeIntegrator<VEC>  dae;

  IndexSet global_partitioning;
  IndexSet partitioning;
  IndexSet relevant_partitioning;

  bool adaptive_refinement;
  bool use_direct_solver;
};

