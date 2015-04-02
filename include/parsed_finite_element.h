#ifndef __dealii_parsed_finite_element_h
#define __dealii_parsed_finite_element_h

#include <deal.II/base/parameter_handler.h>
#include <deal.II/fe/fe.h>

#include "parameter_acceptor.h"

using namespace dealii;

/**
 * Parsed FiniteElement. Read from a parameter file the name of a
 * finite element, generate one for you, and return a pointer to it.
 *
 * The name must be in the form which is returned by the
 * FiniteElement::get_name function, where dimension template
 * parameters <2> etc. can be omitted. Alternatively, the explicit
 * number can be replaced by dim or d. If a number is given, it must
 * match the template parameter of this function.
 *
 * The names of FESystem elements follow the pattern
 * FESystem[FE_Base1^p1-FE_Base2^p2] The powers p1 etc. may either be
 * numbers or can be replaced by dim or d.
 *
 * If no finite element can be reconstructed from this string, an
 * exception of type FETools::ExcInvalidFEName is thrown.
 *
 * The operator() returns a pointer to a newly create finite element. It
 * is in the caller's responsibility to destroy the object pointed to
 * at an appropriate later time.
 *
 * Since the value of the template argument can't be deduced from the
 * (string) argument given to this function, you have to explicitly
 * specify it when you call this function.
 *
 * This function knows about all the standard elements defined in the
 * library. However, it doesn't by default know about elements that
 * you may have defined in your program. To make your own elements
 * known to this function, use the add_fe_name() function.
 */
template <int dim, int spacedim>
class ParsedFiniteElement : public ParameterAcceptor
{
public:
  /**
   * Constructor. Takes a name for the section of the Parameter
   * Handler to use.
   *
   * This class is derived from ParameterAcceptor. Once you
   * constructed an object of this class, if you call
   * ParameterAcceptor::parse_all_parameters(prm), also the
   * parameters of this class will be filled with values from the
   * argument ParameterHandler.
   */
  ParsedFiniteElement (std::string name="");

  /**
   * Declare possible parameters of this class.
   *
   * The only parameter which is declared here is "Finite element
   * space". The parameter must be in the form which is returned by
   * the FiniteElement::get_name function, where dimension template
   * parameters <2> etc. can be omitted. Alternatively, the explicit
   * number can be replaced by dim or d. If a number is given, it
   * must match the template parameter of this function.
   *
   * The names of FESystem elements follow the pattern
   * FESystem[FE_Base1^p1-FE_Base2^p2] The powers p1 etc. may either be
   * numbers or can be replaced by dim or d.
   *
   * If no finite element can be reconstructed from this string, an
   * exception of type FETools::ExcInvalidFEName is thrown.
   *
   * The operator() returns a pointer to a newly create finite element. It
   * is in the caller's responsibility to destroy the object pointed to
   * at an appropriate later time.
   */
  virtual void declare_parameters(ParameterHandler &prm);

  /**
   * Return a pointer to a newly created Finite Element. It will throw
   * an exception if called before any parsing has occured.
   */
  FiniteElement<dim, spacedim> *operator() ();

private:
  /**
   * Finite Element Name.
   */
  std::string fe_name;
};

#endif
