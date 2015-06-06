//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

//
// Simple example of solving the following nonlinear system of
// equations
//
// x(0)^2 + x(1)^2 -1 = 0
//      x(1) - x(0)^2 = 0
//
// using NOX


#include "tests.h"
#include "newton_solver.h"


class SimpleNonlinearProblem:
  public NewtonArgument
{
public:


  unsigned int n_dofs() const
  {
    return 2;
  }

  virtual void output_step(Vector<double> &solution,
                           const unsigned int step_number)
  {
    std::cout<<"Test"<<std::endl;
  }

  /** For newton solver, we need a
  residual function. */
  virtual int residual(Vector<double> &dst,
                       const Vector<double> &src_yy)
  {
    dst(0) = src_yy(0)*src_yy(0) + src_yy(1)*src_yy(1) - 1.0;
    dst(1) = src_yy(1) - src_yy(0)*src_yy(0);
    std::cout<<"In: "<<src_yy<<std::endl;
    std::cout<<"Out: "<<dst<<std::endl;
    return 0;
  }
  /** Jacobian vector product for newton solver. */
  virtual int jacobian(Vector<double> &dst,
                       const Vector<double> &src_yy,
                       const Vector<double> &src)
  {
    double J11 = 2.0*src_yy(0);
    double J12 = 2.0*src_yy(1);
    double J21 = -2.0*src_yy(0);
    double J22 = 1.0;

    dst(0) = J11*src(0) + J12*src(1);
    dst(1) = J21*src(0) + J22*src(1);
    return 0;
  }

  /** Setup Jacobian preconditioner for Newton. */
  virtual int setup_jacobian_prec(const Vector<double> &src_yy)
  {
    return 0;
  }

  /** Jacobian inverse preconditioner
  vector product for newton solver. */
  virtual int jacobian_prec(Vector<double> &dst,
                            const Vector<double> &src_yy,
                            const Vector<double> &src)
  {
    double P11 = 2.0*src_yy(0);
    double P22 = 1.0;

    dst(0) = src(0)/P11;
    dst(1) = src(1)/P22;
    return 0;
  }

  /** Jacobian preconditioner
  vector product for newton solver. */
  virtual int jacobian_prec_prod(Vector<double> &dst,
                                 const Vector<double> &src_yy,
                                 const Vector<double> &src)
  {
    double P11 = 2.0*src_yy(0);
    double P22 = -2.0*src_yy(0);

    dst(0) = src(0)*P11;
    dst(1) = src(1)*P22;
    return 0;
  }

};

int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize init(argc, argv, numbers::invalid_unsigned_int);
  initlog();
  deallog << 0 << std::endl;

  SimpleNonlinearProblem simple_problem;
  NewtonSolver solver(simple_problem);

  Vector<double> solution(2);
  solution(0) = 1.0;
  solution(1) = 1.0;

  solver.solve(solution,30);

  deallog<<"Computed Solution: {"<<solution(0)<<","<<solution(1)<<"}"<<std::endl;
}
