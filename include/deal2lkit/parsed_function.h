#ifndef __dealii_sak_parsed_function_h
#define __dealii_sak_parsed_function_h

#include <deal.II/base/parsed_function.h>
#include <deal2lkit/parameter_acceptor.h>

using namespace dealii;

/**
 * A SAK wrapper for dealii::Functions::ParsedFunction. The template
 * integegers specify the dimension of points this function accepts,
 * and the number of components.
 */
template<int dim, int ncomponents=1>
class ParsedFunction : public ParameterAcceptor, public Functions::ParsedFunction<dim>
{
public:
  /**
   * Constructor: takes an optional name for the section. If the
   * optional expression string is given, than it is used to set the
   * expression as soon as the parameters are declared.
   */
  ParsedFunction(const std::string &name="",
                 const std::string &default_exp="",
                 const std::string &default_const="");

  /**
   * Calls the underlying function of ParsedFunction.
   */
  virtual void declare_parameters(ParameterHandler &prm);

  /**
   * Calls the underlying function of ParsedFunction.
   */
  virtual void parse_parameters(ParameterHandler &prm);


private:
  /**
   * Default expression of this function. "
   */
  const std::string default_exp;
  const std::string default_const;
};

// ============================================================
// Explicit template functions
// ============================================================

template<int dim, int ncomponents>
ParsedFunction<dim, ncomponents>::ParsedFunction(const std::string &name,
                                                 const std::string &default_exp,
                                                 const std::string &default_const) :
  ParameterAcceptor(name),
  Functions::ParsedFunction<dim>(ncomponents),
  default_exp(default_exp),
  default_const(default_const)
{}


template<int dim, int ncomponents>
void ParsedFunction<dim, ncomponents>:: declare_parameters(ParameterHandler &prm)
{
  Functions::ParsedFunction<dim>::declare_parameters(prm, ncomponents);
  if (default_exp != "")
    prm.set("Function expression", default_exp);
  if (default_const != "")
    prm.set("Function constants", default_const);
}


template<int dim, int ncomponents>
void ParsedFunction<dim, ncomponents>:: parse_parameters(ParameterHandler &prm)
{
  Functions::ParsedFunction<dim>::parse_parameters(prm);
}
#endif
