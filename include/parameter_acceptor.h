#ifndef __dealii_parameter_acceptor_h
#define __dealii_parameter_acceptor_h

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/smartpointer.h>
#include <boost/any.hpp>
#include <typeinfo>
#include <string>
#include <iostream>


using namespace dealii;

/**
 * A parameter acceptor base class. This class is used to define a
 * public interface for classes wich use ParameterHandler to obtain
 * class parameters. The only function that derived classes should
 * overload is declare_parameters(). Derived classes are required to
 * use the add_parameter() function inside the declare_paramters()
 * function. If they do so, then
 * ParameterAcceptor::parse_all_parameters() will automatically
 * populate the variables with the parsed parameters. Failure to use
 * the add_parameter() function will result in the user having to call
 * ParameterHandler::get() functions to populate the variables
 * manually.
 *
 * The rationale behind this class, is that we want every class to be
 * able to declare parameters only once, and then automatically
 * populate variables with the values from the parameter file, without
 * having to do so manually.
 */
class ParameterAcceptor : public Subscriptor
{
public:
  /**
   * The constructor adds this class to the list of acceptors. If a
   * section name is specified, then this is used to scope the
   * parameters in the given section, otherwise a pretty printed
   * version of the derived class is used.
   */
  ParameterAcceptor(const std::string section_name="");

  /**
   * The destructor sets to zero the pointer relative to this index,
   * so that it is safe to destroy the mother class.
   */
  ~ParameterAcceptor();


  /**
   * Call declare_all_parameters(), read filename (if it is present as
   * input parameter) and parse_all_parameters() on the static member
   * prm.
   */
  static void initialize(const std::string filename="");

  static void clear();
  /**
   * Parse the parameter file. This function enters the subsection
   * returned by get_section_name() for each derived class, and parse
   * all parameters that were added using add_parameter().
   */
  virtual void parse_parameters(ParameterHandler &prm);

  /**
   * Parse parameter call back. This function is called at the end
   * of parse_parameters, to allow users to process their parameters
   * right after they have been parsed. The default implementation
   * is empty.
   *
   * You can use this function, for example, to create a quadrature
   * rule after you have read how many quadrature points you wanted
   * to use from the parameter file.
   */
  virtual void parse_parameters_call_back();


  /**
   * Generate entries in the given parameter file. Derived classes
   * need to overload this one. If you want to make sure the
   * automatic assignement of variables work, you should fill this
   * function only with calls to the add_parameter() method.
   */
  virtual void declare_parameters(ParameterHandler &prm) = 0;


  /**
   * Parse the given ParameterHandler. This function enters the
   * subsection returned by get_section_name() for each derived class,
   * and parses all parameters that were added using add_parameter().
   */
  static void parse_all_parameters(ParameterHandler &prm);


  /**
   * Initialize the ParameterHandler with all derived classes
   * parameters.This function enters the subsection returned by
   * get_section_name() for each derived class, and declares all
   * parameters that were added using add_parameter().
   */
  static void declare_all_parameters(ParameterHandler &prm);


  /**
   * Return the section name of this class. If a name was provided
   * at construction time, then that name is returned, otherwise it
   * returns the name of this class, pretty printed.
   */
  std::string get_section_name() const;


  /**
   * Add a parameter the given parameter list. A pointer to the
   * parameter is stored, so that every time the default
   * parse_parameters() function is called, this parameter is
   * updated with the value contained in prm.
   */
  template <class T>
  void add_parameter(ParameterHandler &prm, T *parameter,
                     const std::string &entry, const std::string &default_value,
                     const Patterns::PatternBase &pattern=Patterns::Anything(),
                     const std::string &documentation=std::string())
  {
    prm.declare_entry(entry, default_value, pattern, documentation);
    parameters[entry] = boost::any(parameter);
  }

  /**
   * Static parameter. This is used if the user does not provide one.
   */
  static ParameterHandler prm;

private:
  /**
   * A list containing all constructed classes of type
   * ParameterAcceptor.
   */
  static std::vector<SmartPointer<ParameterAcceptor> > class_list;

  /** The index of this specific class within the class list. */
  const unsigned int acceptor_id;

  /**
   * A map of parameters that are initialized in this class with the
   * function add_parameter.
   */
  std::map<std::string, boost::any> parameters;

protected:
  /** The subsection name for this class. */
  const std::string section_name;
};
#endif
