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

#ifndef _d2k_heat_ida_h
#define _d2k_heat_ida_h

#include <deal2lkit/config.h>

#ifdef D2K_WITH_SUNDIALS
#include <deal.II/base/timer.h>
#include <deal.II/base/index_set.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/numerics/error_estimator.h>

#include <deal2lkit/sundials_interface.h>
#include <deal2lkit/ida_interface.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_grid_generator.h>
#include <deal2lkit/parsed_finite_element.h>
#include <deal2lkit/error_handler.h>
#include <deal2lkit/parsed_function.h>
#include <deal2lkit/parsed_data_out.h>
#include <deal2lkit/parsed_dirichlet_bcs.h>
#include <deal2lkit/parsed_solver.h>

#include <fstream>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


using namespace dealii;
using namespace deal2lkit;

typedef TrilinosWrappers::MPI::Vector VEC;
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
                                     const unsigned int step_number,
                                     const double h,
                                     VEC &solution,
                                     VEC &solution_dot);

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
  void refine_mesh ();
  void make_grid_fe();
  void setup_dofs (const bool &first_run=true);

  void assemble_jacobian_matrix (const double t,
                                 const VEC &y,
                                 const VEC &y_dot,
                                 const double alpha);

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

  ErrorHandler<dim>       eh;
  ParsedGridGenerator<dim,dim>   pgg;
  ParsedFiniteElement<dim,dim> fe_builder;

  ParsedFunction<dim, 1>        exact_solution;
  ParsedFunction<dim, 1>        forcing_term;

  ParsedFunction<dim, 1>        initial_solution;
  ParsedFunction<dim, 1>        initial_solution_dot;
  ParsedDirichletBCs<dim,dim,1> dirichlet_bcs;

  ParsedDataOut<dim, dim>       data_out;

  ParsedSolver<VEC>             Ainv;

  IDAInterface<VEC>  dae;

  IndexSet global_partitioning;
  IndexSet partitioning;
  IndexSet relevant_partitioning;

  bool adaptive_refinement;
  bool use_direct_solver;
  bool use_space_adaptivity;
  double kelly_threshold;
  int max_cells;
  double top_fraction;
  double bottom_fraction;
  double diffusivity;
};

#endif

#endif
