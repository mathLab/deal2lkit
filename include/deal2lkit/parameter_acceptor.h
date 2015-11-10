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

/**
@{

@mainpage deal2lkit documentation
@section intro Introduction

\dk is a collection of modules and classes for the general purpose
finite element library \dealii. Its principal aim is to provide a
high level interface, controlled via parameter files, for those
steps that are common in all finite element programs: mesh
generation, selection of the finite element type, application of
boundary conditions and many others. Each module can be used as a
building block independently on the others, and can be integrated
in existing finite element codes based on \dealii, drastically
reducing the size of programs, rendering their use automatically
parametrised, and reducing the overall time-to-market of finite
element programming. Moreover, \dk features interfaces with the
\sundials library (SUite of Nonlinear and DIfferential/ALgebraic
equation Solvers) and with the \assimp library (Open Asset Import
Library).

The \dk library is released under the GNU Lesser General Public
License (LGPL) and can be retrieved from the \dk repository
https://github.com/mathLab/deal2lkit.

@section motivation Motivation

The solution of partial differential equations by means of a finite
element method always requires at least the following steps:

- generation of a geometrical grid to represent the domain of the simulation;
- definition of the discrete functional space for the solution;
- application of proper boundary conditions;
- actual solution of the algebraic problem;
- post-processing of the result (data output and error analysis).

Such a structure usually implies that different problems share a
considerable amount of code. A natural response to such common
background lies in the use of open source libraries as building
blocks for advanced numerical solvers. The general purpose finite
element library \dealii is one of the most successful libraries of
this kind, and allows considerable simplification when writing
complex finite element codes.

The \dealii library has been written with generality in mind, and
allows the solution of several classes of finite element problems. Its
flexibility can be attributed to the granularity and modularity of the
code base, in which only the building blocks of finite element codes
are programmed, and the semantic for the solution of an actual problem
is left to users of the library. This approach has the advantage that
\dealii can be used to solve virtually any problem that can be written
into a partial differential equation, but leaves to the user the
burden to stich together the various building blocks.  A typical
approach is to start from one of the many example programs that the
library comes with (more than 50), and modify it to suite the needs of
the user. While the approach copy-modify-run may be well suited for a
single person working on a single project, it falls short when one
wishes to reuse the same code base to solve possibly very different
problems. The biggest difficulty comes from the fact that most of the
tasks above have slightly different specialisations depending on the
problem at hand. These specialisation are usually difficult to
generalize, since they depend, for example, on the number of variables
of a problem, the types of boundary conditions one would like to
impose, or the type of norm one would like to use when computing
errors during the post-processing phase of a program.

\dk is a **library of modules** built on top of \dealii that
drastically reduces the amount of repeated lines of code between
different projects, by introducing an extensive use of parameter files
into every step of a general finite element code.

\dk features also interfaces for other scientific libraries in order
to tackle problems of increasing difficulties. So far we have
constructed convenience wrappers around the following external
libraries:

- \sundials
- \assimp

\dk is distributed under the free GNU Lesser General Public License
(LGPL) and is available from the \dk repository at
https://github.com/mathLab/deal2lkit. The library is tested by means
of the continuous integration service hosted by Travis CI
(https://travis-ci.org/).

@section overview Modules overview

@subsection parameter_acceptor ParameterAcceptor: the base of all deal2lkit classes

In general, a *parameter file* is used to steer the execution of a
program at run time, without the need to recompile the executable,
with clear advantages in terms of **human-time**.

In the \dealii library, reading and writing parameter files is done
through the ParameterHandler class, that provides a standard
interface to an input file that can be used to feed run-time
parameters to a program, such as time step sizes, geometries, right
hand sides, etc.

\dealii supports the standard `xml` or `JSON` formats, or a
custom text format which resemble bash files with support for
sections, as in the following example:
\code{bash}
subsection Nonlinear solver
  set Nonlinear method = Gradient
  # this is a comment
  subsection Linear solver
    set Solver                       = CG
    set Maximum number of iterations = 30
  end
end
\endcode

Typically, the following four steps are required to let a program use a parameter file:

- make sure that the program knows what entries will be in the file;
- create a parameter file with default values if one does not
  exist;
- parse all entries of the file (possibly raising exceptions if
  the entries were not previously declared, or if the parsed entries
  contain illegal values);
- assign the parsed entries to local variables of the program.

The ParameterHandler class of the \dealii library provides
facilities to perform the above four steps, through the following
methods:

- ParameterHandler::print_parameters()
- ParameterHandler::read_input()
- ParameterHandler::declare_entry()
- ParameterHandler::get()

In large programs, where the number of parameters easily exceeds
hundreds of entries, managing the above four actions for different
classes is far from trivial. The \dealii documentation advocates the
creation of a class that would store all parameters of the problem,
with two methods:

- `declare_parameters(prm)`
- `parse_parameters(prm)`  or `get_parameters(prm)`

that should be called by the program before writing or reading a
parameter file, and right after having read the parameter file into an
object `prm` of type ParameterHandler

Such an approach has the advantage that bookkeeping is simple, if
compared to a scattered approach where each class keeps track of its
own parameters, but it suffers one big draw back: it is not reusable for
problems of different type and it has still the defect that one has to
separate declaration and recovery of each parameter, as in the
following short example:

\code
 void NonLinEq::declare_parameters (ParameterHandler &prm) {
  prm.enter_subsection ("Nonlinear solver");
  {
    prm.declare_entry ("Nonlinear method",
                       "Newton-Raphson",
                       ParameterHandler::RegularExpressions::Anything());
    eq.declare_parameters (prm);
  }
  prm.leave_subsection ();
}
\endcode

The complementary part of this code is contained in the
`parse_parameters` method, which actually fills the values of the
local variables.
\code
void NonLinEq::parse_parameters (ParameterHandler &prm) {
  prm.enter_subsection ("Nonlinear solver");
  std::string method = prm.get ("Nonlinear method");
  eq.parse_parameters (prm);
  prm.leave_subsection ();
}
\endcode

According to the proposed design in the \dealii documentation, such
separation is necessary (with a consequent proliferation of several
places where one has to keep track of what variables have been
declared and what variables have been assigned locally) since the
declaration, reading and writing of a parameter file, and the
assignment to local variables have to be done *exactly* in this
sequence.


\dk implements a **global subscription mechanism** and a **local
subscription mechanism** through the base class ParameterAcceptor,
which maintains compatibility with all classes written following the
\dealii suggested construction, and provides an additional method
which removes the necessity to split the declaration and parsing of
parameters.

The global subscription mechanism is such that whenever a class that
was derived by ParameterAcceptor is constructed, a static registry in
the base class is updated with a pointer to the derived class. Such
registry is traversed upon invocation of the single function
ParameterAcceptor::initialize(file.prm) which in turn calls the method
ParameterAcceptor::declare_parameters() for each of the registered
classes, reads the file `file.prm`, (creating it first with default
values if it does not exist) and subsequently calls the method
ParameterAcceptor::parse_parameters(), again for each of the
registered classes.

@subsection grid_generator ParsedGridGenerator

Of the basic steps for any finite element code, mesh generation and
mesh import are among those tasks which are almost equal in every user
code. ParsedGridGenerator is \dk interface to a collection of \dealii
classed dedicated to creating, reading, and writing a Triangulation to
and from files. ParsedGridGenerator is a wrapper, derived from
ParameterAcceptor, to the following methods and classes:

- GridGenerator: all meshes that \dealii can generate are available
  by selecting their name in a parameter file;

- GridIn: all formats that \dealii can read are available, by
  selecting `file` as the mesh to generate, and then specifying an
  input file name;

- GridOut: selecting a non empty output file name one can create also
  a file containing the Triangulation in any of the output format
  supported by GridOut, by calling the ParsedGridGenerator::write()
  method.
@}
*/


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


D2K_NAMESPACE_CLOSE

#endif

