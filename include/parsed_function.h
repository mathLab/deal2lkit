#ifndef __dealii_sak_parsed_function_h
#define __dealii_sak_parsed_function_h

#include <deal.II/base/parsed_function.h>
#include "parameter_acceptor.h"

using namespace dealii;

/**
 * A SAK wrapper for dealii::Functions::ParsedFunction. The template
 * integegers specify the dimension of points this function accepts,
 * and the number of components. 
 */
template<int dim, int ncomponents=1>
class ParsedFunction : public ParameterAcceptor, Functions::ParsedFunction<dim>
{
public:
  /**
   * Constructor: takes an optional name for the section.
   */ 
  ParsedFunction(std::string name="");

  /**
   * Calls the underlying function of ParsedFunction.
   */
  virtual void declare_parameters(ParameterHandler &prm);
  
};

// ============================================================
// Explicit template functions
// ============================================================

template<int dim, int ncomponents>
ParsedFunction<dim, ncomponents>::ParsedFunction(std::string name) :
  ParameterAcceptor(name),
  Functions::ParsedFunction<dim>(ncomponents)
{}


template<int dim, int ncomponents>
void ParsedFunction<dim, ncomponents>:: declare_parameters(ParameterHandler &prm) {
  Functions::ParsedFunction<dim>::declare_parameters(prm, ncomponents);
}


#endif
