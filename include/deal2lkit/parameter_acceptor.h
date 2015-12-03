//-----------------------------------------------------------
//
//    Copyright (C) 2014 - 2015 by the deal2lkit authors
//
//    This file is part of the deal2lkit library.
//
//    The deal2lkit library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal2lkit distribution.
//
//-----------------------------------------------------------

#ifndef _d2k_parameter_acceptor_h
#define _d2k_parameter_acceptor_h

#include <deal2lkit/config.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <boost/any.hpp>
#include <typeinfo>



using namespace dealii;

D2K_NAMESPACE_OPEN

/**
 * A parameter acceptor base class. This class is used to define a
 * public interface for classes wich want to use a ParameterHandler to
 * handle parameters. Such basic interface provides two subscription
 * mechanisms: a **global subscription mechanism** and a **local
 * subscription mechanism**.
 *
 * The global subscription mechanism is such that whenever a class
 * that was derived by ParameterAcceptor is constructed, a static
 * registry (ParameterAcceptor::class_list) in the base class is
 * updated with a pointer to the derived class. Such registry is
 * traversed upon invocation of the single function
 * ParameterAcceptor::initialize(file.prm) which in turn calls the
 * method ParameterAcceptor::declare_parameters() for each of the
 * registered classes, reads the file `file.prm`, (creating it first
 * with default values if it does not exist) and subsequently calls
 * the method ParameterAcceptor::parse_parameters(), again for each of
 * the registered classes. The method log_info() can be used to
 * extract informations about the classes that have been derived from
 * ParameterAcceptor, and that will be parsed when calling
 * ParameterAcceptor::initialize().
 *
 * ParameterAcceptor conforms to the standard advocated in the \dealii
 * documentation, and it has a pure virtual method
 * ParameterAcceptor::declare_parameters and a virtual method
 * ParameterAcceptor::parse_parameters which can be overloaded as the
 * user whishes. However, the base class also has a default
 * implementation of parse_parameters which exploits a **local
 * subscription mechanism** by storing in a local registry
 * (ParameterAcceptor::parameters) a pointer to all variables that
 * were declared through the ParameterAcceptor::add_parameter
 * method. Such method has the same syntax of the
 * ParameterHandler::add_entry method, with the addition of two
 * arguments: a ParameterHandler object on which
 * ParameterHandler::add_entry will be called, and a reference to the
 * variable that should hold the entry when a ParameterHandler::get_*
 * methods are called. Such variable is stored in
 * ParameterAcceptor::parameters (local to the class instantiation)
 * which is traversed by the default implementation of
 * ParameterAcceptor::parse_parameters. Specialized
 * implementations are provided for the most commonly used variable
 * types.
 *
 * The only function that derived classes should overload is
 * declare_parameters(). Derived classes are required to use the
 * add_parameter() function inside the declare_paramters()
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
 *
 * A typical usage of this chain of classes is the following:
 *
 * @code
 * // This is your own class, derived from ParameterAcceptor
 * class MyClass : public ParameterAcceptor {
 * virtual declare_parameters(ParameterHandler &prm) {
 *  add_parameters(prm, &member_var, "A param", "Default value");
 * }
 *
 * ...
 *
 * int main() {
 *  // Make sure you build your class BEFORE calling
 *  // ParameterAcceptor::initialize()
 *  MyClass class;
 *
 *  // With this call, all derived classes will have their
 *  // parameters initialized
 *  ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 *
 * As a convention, in \dk all modules derived from ParameterAcceptor
 * implement a default constructor which takes one or more optional
 * arguments. The first optional argument is always the name of the
 * section in the parameter file that the derived class should use to
 * fill its local variables. By default, the utility function
 * deal2lkit::type is used to fill the section name with a human
 * readable version of the class name itself. For example, according
 * to the above example, there will be a section named
 * MyClass. Next options are the default values that will be
 * written in the parameter file. There is no need to specify these
 * options, as the user can always change the content of the file to
 * make sure that the right parameters are used at run time, but this
 * possibility allows one to design a program that does something
 * sensible on the first run, without having to change any parameter
 * file.
 */
class ParameterAcceptor : public Subscriptor
{
public:
  /**
   * The constructor adds derived classes to the list of acceptors. If
   * a section name is specified, then this is used to scope the
   * parameters in the given section, otherwise a pretty printed
   * version of the derived class is used.
   */
  ParameterAcceptor(const std::string section_name="");

  /**
   * The destructor sets to zero the pointer relative to this index,
   * so that it is safe to destroy the mother class.
   */
  virtual ~ParameterAcceptor();


  /**
   * Call declare_all_parameters(), read filename (if it is present as
   * input parameter) and parse_all_parameters() on the static member
   * prm. If outfilename is present, then write the content that was
   * read in to the outfilename. This is useful to get all used
   * parameters, and not only those that were set in the input file.
   */
  static void initialize(const std::string filename="",
                         const std::string outfilename="");

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
  static void parse_all_parameters(ParameterHandler &prm=ParameterAcceptor::prm);


  /**
   * Print information about all stored classes.
   */
  static void log_info();


  /**
   * Initialize the ParameterHandler with all derived classes
   * parameters.This function enters the subsection returned by
   * get_section_name() for each derived class, and declares all
   * parameters that were added using add_parameter().
   */
  static void declare_all_parameters(ParameterHandler &prm=ParameterAcceptor::prm);

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
                     const std::string &entry,
                     const std::string &default_value,
                     const Patterns::PatternBase &pattern=Patterns::Anything(),
                     const std::string &documentation=std::string())
  {
    AssertThrow(std::is_const<T>::value == false,
                ExcMessage("You tried to add a parameter using a const "
                           "variable. This is not allowed, since these "
                           "variables will be filled later on when "
                           "parsing the parameter."));

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


D2K_NAMESPACE_CLOSE

#endif

