#ifndef __dealii_parameter_acceptor_h
#define __dealii_parameter_acceptor_h

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/smartpointer.h>
#include <typeinfo>

using namespace dealii;

/** 
  A parameter acceptor base class. This is a pure virtual class, used
  to define a public interface for classes with declare_parameters()
  and parse_parameters() methods.
*/
class ParameterAcceptor : public Subscriptor
{
public:
  /** The constructor adds this class to the list of acceptors. */
  ParameterAcceptor(const std::string section_name);

  /** The destructor sets to zero the pointer relative to this index,
      so that it is safe to destroy the mother class. */
  ~ParameterAcceptor();
  
  /** Parse the parameter file. */
  virtual void parse_parameters(ParameterHandler &prm) = 0;

  /** Generate entries in the given parameter file. */
  virtual void declare_parameters(ParameterHandler &prm) = 0;

  /** Call each of the parse_parameters in the class_list vector. */
  void parse_all_parameters(ParameterHandler &prm);
  
  /** Call each of the declare_parameters in the class_list vector. */
  void declare_all_parameters(ParameterHandler &prm);

private: 
  /** A list containing all constructed classes of type
      ParameterAcceptor. */
  static std::vector<SmartPointer<ParameterAcceptor> > class_list;

  /** The index of this specific class within the class list.*/
  const unsigned int acceptor_id;    

protected:
  /** The subsection name for this class. */
  const std::string section_name;
};
#endif
