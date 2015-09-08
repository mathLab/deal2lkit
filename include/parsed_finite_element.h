#ifndef __dealii_parsed_finite_element_h
#define __dealii_parsed_finite_element_h

#include <deal.II/base/parameter_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>

#include "parameter_acceptor.h"
#include "utilities.h"

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
template <int dim, int spacedim=dim>
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
   *
   * The optional parameters specify the FiniteElement name, the
   * component names, the allowed number of components and the system
   * and preconditioner couplings. This class will throw an exception
   * if the number of components is set to something different than
   * zero and the corresponding FiniteElement does not match this
   * number of components. If n_components is left to 0 (the default
   * value), then any FiniteElement can be generated, with arbitrary
   * numbers of components.
   */
  ParsedFiniteElement (const std::string &name="",
                       const std::string &default_fe="FE_Q(1)",
                       const std::string &default_component_names="u",
                       const unsigned int n_components=0,
                       const std::string &default_coupling="",
                       const std::string &default_preconditioner_coupling="");

  /**
   * Declare possible parameters of this class.
   *
   * The first parameter which is declared here is "Finite element
   * space". The parameter must be in the form which is returned by
   * the FiniteElement::get_name function, where dimension template
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
   * The second parameter is "Block names". This is comma separeted
   * list of component names which identifies the Finite Elemenet. If
   * a name is repeated, then that component is assumed to be part of
   * a vector field or Tensor field, and is treated as a single
   * block. User classes can use this information to construct block
   * matrices and vectors, or to group solution names according to
   * components. For example, a Stokes problem may have "u,u,p" for
   * dim = 2 or "u, u, u, p" for dim = 3.
   *
   * The third parameter defined by this function is "Block coupling",
   * which is semicolumn separated list of comma separated indices
   * from 0 to 2 (included), representing a table of coupling between
   * the components (or the blocks). The default for this parameter is
   * the empty string, which is interpreted as "all components couple
   * with all other components". You can specify a block wise coupling
   * (with two blocks, this would be a 2x2 matrix) or a component wise
   * coupling (if a block is made of more than one component, but not
   * all components of the same block couple with the other
   * components). The numbers are interpreted as: 0=DoFTools::none,
   * 1=DoFTools::always, 2=DoFTools::nonzero.
   */
  virtual void declare_parameters(ParameterHandler &prm);

  /**
   * Return a shared pointer to a newly created Finite Element. It
   * will throw an exception if called before any parsing has occured,
   * or if the number of components of the generated finite element is
   * different from the number of components given at construction
   * time.
   */
  FiniteElement<dim, spacedim> *operator() () const;

  /**
   * Fill information about blocks after parsing the parameters.
   */
  virtual void parse_parameters_call_back();

  /**
   * Return the component names for this Finite Element.
   */
  std::string get_component_names() const;

  /**
   * Return the blocking of the components for this finite
   * element. This is what's needed by the block renumbering
   * algorithm.
   */
  std::vector<unsigned int> get_component_blocks() const;


  /**
   * Return the block names for this Finite Element. This is the same
   * as std::unique(get_component_names().begin(),
   * get_component_names().end())
   */
  std::string get_block_names() const;

  /**
   * Return the number of components of the Finite Element.
   */
  unsigned int n_components() const;

  /**
   * Return the number of blocks of the Finite Element, i.e.,
   * the number of variables. For example, simple Heat equation
   * has 1 block, Navier-Stokes 2 blocks (u and p).
   */
  unsigned int n_blocks() const;

  /**
   * Return the coupling of the Finite Element for the system matrix.
   */
  const Table<2,DoFTools::Coupling> &get_coupling() const;


  /**
   * Return the coupling of the Finite Element for the preconditioner.
   */
  const Table<2,DoFTools::Coupling> &get_preconditioner_coupling() const;


private:
  /**
   * Convert integer table into a coupling table.
   */
  Table<2, DoFTools::Coupling> to_coupling(const std::vector<std::vector<unsigned int> > &table) const;

  /**
   * Number of components of this FiniteElement. If you want to allow
   * for arbitrary components, leave this to its default value 0.
   */
  const unsigned int _n_components;

  /**
   * Finite Element Name.
   */
  std::string fe_name;

  /**
   * Default component names.
   */
  std::string default_component_names;

  /**
   * Default system coupling.
   */
  std::string default_coupling;

  /**
   * Default preconditioner coupling.
   */
  std::string default_preconditioner_coupling;

  /**
   * Block names. This is comma separeted list of component names
   * which identifies the Finite Elemenet. If a name is repeated, then
   * that component is assumed to be part of a vector field or Tensor
   * field, and is treated as a single block. User classes can use
   * this information to construct block matrices and vectors, or to
   * group solution names according to components. For example, a
   * Stokes problem may have "u,u,p" for dim = 2 or "u, u, u, p" for
   * dim = 3.
   */
  std::vector<std::string> component_names;

  /**
   * The subdivision, in terms of component indices. This is
   * automatically computed from the the component names.
   */
  std::vector<unsigned int> component_blocks;


  /**
   * The subdivision, in terms of block names. This is automatically
   * computed from the the component names.
   */
  std::vector<std::string> block_names;

  /**
   * System coupling information. This is a readable version of the
   * coupling information. This is transformed after reading the
   * parameter to a table of DoFTools::coupling(), which is needed by
   * the sparsity_pattern builders.
   */
  std::vector<std::vector<unsigned int> > coupling_int;


  /**
   * Preconditioner coupling information. This is a readable version
   * of the coupling information. This is transformed after reading
   * the parameter to a table of DoFTools::coupling(), which is needed
   * by the sparsity_pattern builders.
   */
  std::vector<std::vector<unsigned int> > preconditioner_coupling_int;

  /**
   * Coupling information for the system matrix in a format compatible
   * for sparsity construction.
   */
  Table<2,DoFTools::Coupling> coupling;

  /**
   * Coupling information for the preconditioner in a format
   * compatible for sparsity construction.
   */
  Table<2,DoFTools::Coupling> preconditioner_coupling;
};

#endif
