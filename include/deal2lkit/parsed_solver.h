#ifndef __dealii_sak_parsed_solver_h
#define __dealii_sak_parsed_solver_h

#include <deal2lkit/config.h>
#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CXX11

#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_qmrs.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/linear_operator.h>

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/utilities.h>

using namespace dealii;

template <typename VECTOR>
std::function<void(VECTOR &, bool)> default_reinit()
{
  return [](VECTOR &, bool)
  {
    Assert(false, ExcInternalError("It seems you really need to set a ReinitFunction for your "
                                   "operator. Try avoiding PackagedOperator if you don't want "
                                   "to implement this function."));
  };
}

/**
 * A solver selector which uses parameter files to choose between
 * different options. This object is a LinearOperator which can be
 * called in place of the inverse of another LinearOperator.
 *
 * Example usage is the following:
 *
 * @code
 * ParsedSolver<VEC> Ainv;
 * ParameterAcceptor::initialize(...);
 *
 * Ainv.op   = linear_operator<VEC>(A);
 * Ainv.prec = linear_operator<VEC>(sompreconditioners);
 *
 * x = Ainv*b;
 *
 * @endcode
 */
template<typename VECTOR>
class ParsedSolver : public LinearOperator<VECTOR, VECTOR>,  public ParameterAcceptor
{
public:
  /**
   * Constructor. Build the inverse of an Operator using a parameter
   * file. A section name can be specified, the solver type, the
   * maximum number of iterations, and the reduction to reach
   * convergence. If you know in advance the operators this object
   * will need, you can also supply them here. They default to the
   * identity, and you can assign them later by setting op and prec.
   */
  ParsedSolver(const std::string &name="",
               const std::string &default_solver="cg",
               const long int iter=1000,
               const double reduction=1e-8,
               const LinearOperator<VECTOR> &op=
                 identity_operator<VECTOR>(default_reinit<VECTOR>()),
               const LinearOperator<VECTOR> &prec=
                 identity_operator<VECTOR>(default_reinit<VECTOR>()));

  /**
   * Declare solver type and solver options.
   */
  virtual void declare_parameters(ParameterHandler &prm);

  /**
   * Parse solver type and solver options.
   */
  virtual void parse_parameters(ParameterHandler &prm);

  /**
   * Initialize internal variables.
   */
  virtual void parse_parameters_call_back();


  /**
   * The Operator this solver use. You can assign a new one at
   * construction time, or by this->op = some_new_op. By default it is
   * the identity operator.
   */
  LinearOperator<VECTOR> op;

  /**
   * The preconditioner used by this solver. You can assign a new one
   * at construction time, or by this->op = some_new_op.  By default
   * it is the identity operator.
   */
  LinearOperator<VECTOR> prec;

  /**
   * ReductionControl. Used internally by the solver.
   */
  ReductionControl control;
private:
  /**
   * Store a shared pointer, and intilize the inverse operator.
   */
  template<typename MySolver >
  void initialize_solver(MySolver *);

  /**
   * Solver name."
   */
  std::string solver_name;

  /**
   * Default number of maximum iterations required to succesfully
   * complete a solution step.
   */
  long int max_iterations;

  /**
   * Default reduction required to succesfully complete a solution
   * step.
   */
  double reduction;

  /**
   * The actual solver.
   */
  shared_ptr<Solver<VECTOR> > solver;
};

// ============================================================
// Explicit template functions
// ============================================================

template<typename VECTOR>
ParsedSolver<VECTOR>::ParsedSolver(const std::string &name,
                                   const std::string &default_solver,
                                   const long int default_iter,
                                   const double default_reduction,
                                   const LinearOperator<VECTOR> &op,
                                   const LinearOperator<VECTOR> &prec) :
  ParameterAcceptor(name),
  op(op),
  prec(prec),
  solver_name(default_solver),
  max_iterations(default_iter),
  reduction(default_reduction)
{}


template<typename VECTOR>
void ParsedSolver<VECTOR>::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &solver_name, "Solver name", solver_name,
                Patterns::Selection("cg|bicgstab|gmres|fgmres|"
                                    "minres|qmrs|richardson"),
                "Name of the solver to use.");

  ReductionControl::declare_parameters(prm);

  prm.set("Max steps", max_iterations);
  prm.set("Reduction", reduction);
}


template<typename VECTOR>
void ParsedSolver<VECTOR>::parse_parameters(ParameterHandler &prm)
{
  ParameterAcceptor::parse_parameters(prm);
  control.parse_parameters(prm);
}


template<typename VECTOR>
template<typename MySolver >
void ParsedSolver<VECTOR>::initialize_solver(MySolver *s)
{
  solver = SP(s);
  (LinearOperator<VECTOR,VECTOR> &)(*this) = inverse_operator(op, *s, prec);
}

template<typename VECTOR>
void ParsedSolver<VECTOR>::parse_parameters_call_back()
{
  if (solver_name == "cg")
    {
      initialize_solver(new SolverCG<VECTOR>(control));
    }
  else if (solver_name == "bicgstab")
    {
      initialize_solver(new SolverBicgstab<VECTOR>(control));
    }
  else if (solver_name == "gmres")
    {
      initialize_solver(new SolverGMRES<VECTOR>(control));
    }
  else if (solver_name == "fgmres")
    {
      initialize_solver(new SolverFGMRES<VECTOR>(control));
    }
  else if (solver_name == "minres")
    {
      initialize_solver(new SolverMinRes<VECTOR>(control));
    }
  else if (solver_name == "qmrs")
    {
      initialize_solver(new SolverQMRS<VECTOR>(control));
    }
  else if (solver_name == "richardson")
    {
      initialize_solver(new SolverRichardson<VECTOR>(control));
    }
  else
    {
      Assert(false, ExcInternalError("Solver should not be unknonw."));
    }
}

#endif // CXX11
#endif
