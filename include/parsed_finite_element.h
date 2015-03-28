#ifndef __dealii_parsed_finite_element_h
#define __dealii_parsed_finite_element_h

#include <deal.II/base/parameter_handler.h>
#include <deal.II/fe/fe.h>

#include "parameter_acceptor.h"

using namespace dealii;

/**
 * Parsed FiniteElement. Read from a parameter file the name of a
 * finite element, and generate a pointer to it.
 */
template <int dim, int spacedim>
class ParsedFiniteElement : public ParameterAcceptor
{
public:
    /**
     * Constructor. Takes a name for the section of the Parameter Handler to use.
     */
    ParsedFiniteElement (std::string name="");

    /**
     * Declare possible parameters of this class.
     */
    virtual void declare_parameters(ParameterHandler &prm);

    /**
     * Return a pointer to a newly created Finite Element. It will throw
     * an exception if called before any parsing has occured.
     */
    FiniteElement<dim, spacedim> * operator() ();

private:
    /**
     * Finite Element Name.
     */
    std::string fe_name;
};

#endif
