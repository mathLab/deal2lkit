//-----------------------------------------------------------
//
//    Copyright (C) 2014 - 2016 by the deal2lkit authors
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
#include <deal2lkit/utilities.h>
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
 * registered classes, reads the file `file.prm` (creating it first
 * with default values if it does not exist) and subsequently calls
 * the method ParameterAcceptor::parse_parameters(), again for each of
 * the registered classes. The method log_info() can be used to
 * extract informations about the classes that have been derived from
 * ParameterAcceptor, and that will be parsed when calling
 * ParameterAcceptor::initialize().
 *
 * ParameterAcceptor conforms to the standard advocated in the \dealii
 * documentation, and it has a virtual method
 * ParameterAcceptor::declare_parameters and a virtual method
 * ParameterAcceptor::parse_parameters which can be overloaded as the
 * user whishes. However, the base class also has a default
 * implementation of parse_parameters which exploits a **local
 * subscription mechanism** by storing in a local registry
 * (ParameterAcceptor::parameters) a pointer to all variables that
 * were declared through the ParameterAcceptor::add_parameter
 * method. Such method has two flavours. In the first we use
 * the same syntax of the
 * ParameterHandler::add_entry method, with the addition of two
 * arguments: a ParameterHandler object on which
 * ParameterHandler::add_entry will be called, and a reference to the
 * variable that should hold the entry when a ParameterHandler::get_*
 * methods are called. The second method instead only takes a reference
 * to the variable and to the text entry in the parameter file (with
 * optional documentation).
 *
 * The variable is stored internally in the
 * ParameterAcceptor::parameters (local to the class instantiation)
 * which is traversed by the default implementation of
 * ParameterAcceptor::parse_parameters. Specialized
 * implementations are provided for the most commonly used variable
 * types.
 *
 * Derived classes are required to use one of the
 * add_parameter() functions, either inside a declare_paramters()
 * function, or in the constructor of the class.
 *
 * In either way, ParameterAcceptor::parse_all_parameters() will automatically
 * populate the variables with the parsed parameters. Failure to use
 * the add_parameter() function will result in the user having to call
 * ParameterHandler::get() functions to populate the variables
 * manually. If some post processing is required on the parsed values,
 * the virtual function ParameterAcceptor::parse_parameters_call_back()
 * can be overridden, which is called just after the parse_parameters()
 * function of each class.
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
 *
 * // The constructor of ParameterAcceptor requires a std::string,
 * // which defines the section name where the parameters of MyClass
 * // will be stored.
 *
 * MyClass(std::string name) :
 *   ParameterAcceptor(name)
 * {}
 *
 * virtual declare_parameters(ParameterHandler &prm) {
 *  add_parameters(prm, &member_var, "A param", "Default value");
 * }
 *
 * ...
 * };
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
 * An even simpler implementation (not conforming to \dealii advocated
 * standards) is given by the following example:
 *
 * @code
 * // Again your own class, derived from ParameterAcceptor
 * class MyClass : public ParameterAcceptor {
 *
 * // We now fill the parameters inside the constructor
 *
 * MyClass(std::string name) :
 *   ParameterAcceptor(name)
 * {
 *  add_parameter(&member_var, "A param", "Documentation");
 *  add_parameter(&another_member_var, "Another param");
 * }
 * ...
 * };
 *
 * int main() {
 *  // Make sure you build your class BEFORE calling
 *  // ParameterAcceptor::initialize()
 *  MyClass class;
 *  ParameterAcceptor::initialize("file.prm");
 *  class.run();
 * }
 * @endcode
 *
 *
 * Parameter files can be organised into section/subsection/subsubsection.
 * To do so, the std::string passed to ParameterAcceptor within the
 * constructor of the derived class needs to contain the separator "/".
 * In fact, "first/second/third/My Class" will organize the parameters
 * as follows
 *
 * @code
 * subsection first
 *   subsection second
 *     subsection third
 *       subsection My Class
 *        ... # all the parameters
 *       end
 *     end
 *   end
 * end
 * @endcode
 *
 * Let's now discuss some cases with increasing complexities in order
 * to understand the best way to manage them with ParameterAcceptor.
 *
 * MyClass is derived from ParameterAcceptor and has a
 * member object that is derived itself from ParameterAcceptor.
 * @code
 * class MyClass : public ParameterAcceptor
 * {
 *   MyClass (std::string name);
 *   virtual void declare_parameters(ParameterHandler &prm);
 * private:
 *   ParsedFunction<dim> function;
 *  ...
 * };
 *
 * MyClass::MyClass(std::string name)
 *  :
 * ParameterAcceptor(name),
 * function("Forcing term")
 * {}
 *
 * void MyClass::declare_parmeters(ParameterHandler &prm)
 * {
 *  // many add_parameter(...);
 * }
 *
 * ...
 *
 * int main()
 * {
 * MyClass mc("My Class");
 *
 * ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 *
 * In this case, the structure of the parameters will be
 * @code
 * subsection Forcing term
 * ... #parameters of ParsedFunction
 * end
 * subsection My class
 * ... #all the parameters of MyClass defined in declare_parameters
 * end
 * @endcode
 * Note that the sections are alphabetically sorted.
 *
 * Now suppose that in the main file we need two or more objects of MyClass
 * @code
 * int main()
 * {
 *  MyClass ca("Class A");
 *  MyClass cb("Class B");
 *  ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 *
 * What we will read in the parameter file looks like
 * @code
 * subsection Class A
 * ...
 * end
 * subsection Class B
 * ...
 * end
 * subsection Forcing term
 * ...
 * end
 * @endcode
 * Note that there is only one section "Forcing term", this is because
 * both objects have defined the same name for the section of their
 * ParsedFunction. There are two strategies to manage this issue. The
 * first one (not recommended) would be to change the name of the section
 * of the ParsedFunction such that it contains also the string passed to
 * the constructor of MyClass:
 * @code
 * MyClass::MyClass(std::string name)
 *  :
 * ParameterAcceptor(name),
 * function(name+" --- forcing term")
 * {}
 * @endcode
 *
 * The other way to proceed (recommended) is to use exploit the /section/subsection
 * approach **in the main class**.
 * @code
 * int main()
 * {
 *  MyClass ca("/Class A/Class");
 *  MyClass cb("/Class B/Class");
 *  ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 * Now, in the parameter file we can find
 * @code
 * subsection Class A
 *   subsection Class
 *   ...
 *   end
 *   subsection Forcing term
 *   ...
 *   end
 * end
 * subsection Class B
 *   subsection Class
 *   ...
 *   end
 *   subsection Forcing term
 *   ...
 *   end
 * end
 * @endcode
 * Note the "/" at the begin of the string name. This is interpreted by
 * ParameterAcceptor like the root folder in Unix systems. This means
 * that the sections "Class A" and "Class B" will not be nested under any
 * section. On the other hand, if the string does not begin with a "/"
 * as in the previous cases (and for the ParsedFunction also in this last
 * example) the section will be created **under the current path**, which
 * depends on the previously defined sections/subsections/subsubsections.
 * Indeed, the section "Forcing term" is nested under "Class A" or "Class B".
 * To make things more clear. let's consider the following two examples
 * @code
 * int main()
 * {
 *  MyClass ca("/Class A/Class");
 *  MyClass cb("Class B/Class");
 *  ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 * The parameter file will have the following structure
 * @code
 * subsection Class A
 *   subsection Class
 *   ...
 *   end
 *   subsection Forcing term
 *   ...
 *   end
 *   subsection Class B
 *     subsection Class
 *     ...
 *     end
 *     subsection Forcing term
 *     ...
 *     end
 *   end
 * end
 * @endcode
 *
 * If instead one of the paths ends with "/" instead of just
 * a name of the class, subsequen classes will be declared
 * under the full path, as if the class name should be interpreted
 * as a directory:
 * @code
 * int main()
 * {
 *  MyClass ca("/Class A/Class/");
 *  MyClass cb("Class B/Class");
 *  ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 * The parameter file will have the following structure
 * @code
 * subsection Class A
 *   subsection Class
 *   ...
 *      subsection Forcing term
 *      ...
 *      end
 *      subsection Class B
 *          subsection Class
 *          ...
 *          end
 *          subsection Forcing term
 *          ...
 *          end
 *      end
 *   end
 * end
 * @endcode
 *
 * As a final remark, in order to allow a proper management of all the
 * sections/subsections, the instantiation of objects and the call to
 * ParameterAcceptor::initialize() **cannot be done in multithread**.
 *
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

  /**
   * Clear class list and global parameter file.
   */
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
   * may want to overload this function. If you want to make sure the
   * automatic assignement of variables work, you should fill this
   * function only with calls to the add_parameter() method.
   *
   * Alternatively, you could insert your calls to the add_parameter()
   * function that takes three arguments directly in the constructor.
   *
   * In this case, the ParameterAcceptor::prm parameter handler is used
   * by default, and you don't need to overload this function. The default
   * implementation in fact does nothing.
   *
   * In general this approach is here to guarantee backward compatibility
   * with the strategy advocated by the \dealii library of splitting declaration
   * and parsing of parameters into two functions.
   */
  virtual void declare_parameters(ParameterHandler &) {};


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
   * Travers all registered classes, and figure out what
   * subsections we need to enter.
   */
  std::vector<std::string> get_section_path() const;

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
    patterns[entry] = SP(pattern.clone());
  }


  /**
   * Add a parameter to the global parameter handler ParameterAcceptor::prm.
   * A pointer to the parameter is stored, so that every time the default
   * parse_parameters() function is called, this parameter is
   * updated with the value contained in prm. The default value of the
   * parameter is taken by the current content of the parameter, and the
   * default Pattern is constructed using the method to_pattern().
   *
   * Notice that this function has a slightly different behaviour with respect
   * to the other add_parameter() method, since it assumes that the global
   * parameter is always in its root section, and therefore before calling
   * prm.add_entry() it will enter in the sections specified by this class,
   * and leave all entered sections after having declared the variable, leaving
   * the parameter in the same state as before, but having inserted the entry
   * in the nested sections returned by get_section_path().
   */
  template <class T>
  void add_parameter(T &parameter,
                     const std::string &entry,
                     const Patterns::PatternBase &pattern=*to_pattern(T()),
                     const std::string &documentation=std::string(),
                     ParameterHandler &prm=ParameterAcceptor::prm)
  {
    AssertThrow(std::is_const<T>::value == false,
                ExcMessage("You tried to add a parameter using a const "
                           "variable. This is not allowed, since these "
                           "variables will be filled later on when "
                           "parsing the parameter."));

    enter_my_subsection(prm);
    prm.declare_entry(entry, to_string(parameter),
                      pattern,
                      documentation);
    leave_my_subsection(prm);
    parameters[entry] = boost::any(&parameter);
    patterns[entry] = SP(pattern.clone());
  }

  /**
   * Make sure we enter the right subsection of the global parameter file.
   * This function should be called when prm is in its root subsection.
   */
  void enter_my_subsection(ParameterHandler &prm);

  /**
   * This function undoes what the enter_my_subsection() function did. It only
   * makes sense if enter_my_subsection() is called before this one.
   */
  void leave_my_subsection(ParameterHandler &prm);

  /**
   * Given a class T, construct its default pattern to be used when declaring
   * parameters.
   */
  template <class T>
  static std::shared_ptr<Patterns::PatternBase> to_pattern(const T &);

  /**
   * Given a string, fill the value of the given parameter.
   */
  template <class T>
  static T to_type(const std::string &);

  /**
   * Given a parameter, return a string containing the given parameter.
   */
  template <class T>
  static std::string to_string(const T &);

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
   * functions add_parameters.
   */
  mutable std::map<std::string, boost::any> parameters;

  /**
   * A map of patterns that are initialized in this class with the
   * functions add_parameters.
   */
  mutable std::map<std::string, std::shared_ptr<Patterns::PatternBase> > patterns;


  /**
   * Separator between section and subsection.
   */
  static const char sep = '/';

protected:
  /** The subsection name for this class. */
  const std::string section_name;
};


D2K_NAMESPACE_CLOSE

#endif

