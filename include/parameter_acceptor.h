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

  The rationale behind this class, is that we want to make sure all
  classes declare and parse their own parameters. In order to ensure
  this, the base class keeps track of pointers to derived classes, and
  two functions are available which collectively call each declare and
  parse parameter function.
*/
class ParameterAcceptor : public Subscriptor
{
public:
  /** The constructor adds this class to the list of acceptors. If a
      section name is specified, then this is used to scope the
      parameters in the given section, otherwise a pretty printed
      version of the derived class is used. */
  ParameterAcceptor(const std::string section_name="");

  /** The destructor sets to zero the pointer relative to this index,
      so that it is safe to destroy the mother class. */
  ~ParameterAcceptor();
  
  /** Parse the parameter file. Derived classes need to overload
      this one. */
  virtual void parse_parameters(ParameterHandler &prm) = 0;

  /** Generate entries in the given parameter file. Derived classes
      need to overload this one. */
  virtual void declare_parameters(ParameterHandler &prm) = 0;

  /** Call each of the parse_parameters in the class_list vector. */
  void parse_all_parameters(ParameterHandler &prm);
  
  /** Call each of the declare_parameters in the class_list vector. */
  void declare_all_parameters(ParameterHandler &prm);

  /** Return the section name of this class. If a name was provided at
      construction time, then that name is returned, otherwise it
      returns the name of this class, pretty printed. */
  std::string get_section_name() const;

private: 
  /** A list containing all constructed classes of type
      ParameterAcceptor. */
  static std::vector<SmartPointer<ParameterAcceptor> > class_list;
  static std::vector<std::string> section_names;

  /** The index of this specific class within the class list.*/
  const unsigned int acceptor_id;    

protected:
  /** The subsection name for this class. */
  const std::string section_name;
};
#endif
