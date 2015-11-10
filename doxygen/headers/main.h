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
