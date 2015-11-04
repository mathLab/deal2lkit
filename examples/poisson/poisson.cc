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

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_grid_generator.h>
#include <deal2lkit/parsed_finite_element.h>
#include <deal2lkit/parsed_data_out.h>
#include <deal2lkit/parsed_dirichlet_bcs.h>
#include <deal2lkit/error_handler.h>
#include <deal2lkit/parsed_solver.h>
#include <deal2lkit/utilities.h>


#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>


#include <fstream>

using namespace deal2lkit;
using namespace dealii;

/** Solve the problem
 *
 *  \f[
 *  -\nabla\cdot(kappa \nabla u) = f
 *  \f]
 */

class PoissonParameters : public ParameterAcceptor
{
public:
  virtual void declare_parameters(ParameterHandler &prm)
  {
    add_parameter(prm, &n_cycles, "Number of cycles", "4");
    add_parameter(prm, &initial_refinement, "Initial refinement", "3");
  }
  unsigned int n_cycles;
  unsigned int initial_refinement;
};

int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  // Change this if you want a three dimensional simulation
  const unsigned int dim = 2;
  const unsigned int spacedim = 2;

  PoissonParameters par;
  ParsedGridGenerator<dim,spacedim> pgg;

  ParsedFiniteElement<dim,spacedim> pfe;
  ParsedDirichletBCs<dim,dim,1> bcs;
  ParsedFunction<spacedim> kappa("Kappa", "1.0");
  ParsedFunction<spacedim> force("Forcing term");
  ParsedFunction<spacedim> exact("Exact solution");


  SparsityPattern sparsity;
  SparseMatrix<double> matrix;
  ConstraintMatrix constraints;
  Vector<double> solution;
  Vector<double> rhs;

  PreconditionIdentity prec;

  ParsedSolver<Vector<double> > inverse("Solver",
                                        "cg",
                                        1000,
                                        1e-8,
                                        linear_operator<Vector<double> >(matrix),
                                        linear_operator<Vector<double> >(matrix, prec));



  shared_ptr<Triangulation<dim,spacedim> > tria;
  shared_ptr<FiniteElement<dim,spacedim> > fe;

  ErrorHandler<1> eh;
  ParsedDataOut<dim,spacedim> pdo;


  // Generate Triangulation, DoFHandler, and FE
  ParameterAcceptor::initialize("poisson.prm", "used_parameters.prm");


  tria = SP(pgg.serial());
  pgg.write(*tria);

  fe = SP(pfe());
  DoFHandler<dim,spacedim> dh(*tria);

  tria->refine_global(par.initial_refinement);
  QGauss<dim> quad(2*fe->degree+1);

  for (unsigned int i=0; i<par.n_cycles; ++i)
    {
      dh.distribute_dofs(*fe);
      std::cout << "Cycle " << i
                << ", cells: " << tria->n_active_cells()
                << ", dofs: " << dh.n_dofs() << std::endl;

      DynamicSparsityPattern dsp(dh.n_dofs());
      DoFTools::make_sparsity_pattern (dh, dsp);
      sparsity.copy_from(dsp);

      matrix.reinit (sparsity);
      solution.reinit (dh.n_dofs());
      rhs.reinit (dh.n_dofs());

      constraints.clear();
      bcs.interpolate_boundary_values(dh, constraints);
      constraints.close();

      MatrixCreator::create_laplace_matrix(dh, quad, matrix, force, rhs, &kappa);

      constraints.condense(matrix, rhs);


      // solution = inverse*rhs;

      SparseDirectUMFPACK solver;
      solver.initialize(matrix);
      solver.vmult(solution, rhs);
      constraints.distribute(solution);

      pdo.prepare_data_output(dh, std::to_string(i));
      pdo.add_data_vector(solution, "solution");
      pdo.write_data_and_clear();

      eh.error_from_exact(dh, solution, exact);

      tria->refine_global(1);
    }

  eh.output_table();
  return 0;
}
