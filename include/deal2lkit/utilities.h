//-----------------------------------------------------------
//
//    Copyright (C) 2015 by the deal2lkit authors
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

#ifndef _d2k_utilities_h
#define _d2k_utilities_h

#include <deal2lkit/config.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/timer.h>
#include <deal.II/fe/fe.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>

#include <typeinfo>
#include <cxxabi.h>
#include <sstream>
#include <sys/ioctl.h>    // to know the number of cols and rows of a shell
#include <chrono>         // for TimeUtilities std::chrono
#include <stdio.h>

#include <deal.II/lac/block_vector.h>
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#ifdef DEAL_II_WITH_TRILINOS
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <Teuchos_TimeMonitor.hpp>
#endif
#ifdef DEAL_II_WITH_PETSC
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#endif
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#include <deal.II/base/utilities.h>

#ifdef D2K_WITH_SUNDIALS
#include <nvector/nvector_serial.h>
#ifdef DEAL_II_WITH_MPI
#include <nvector/nvector_parallel.h>
#endif
#endif


#include <deal.II/base/index_set.h>
using namespace dealii;
using std_cxx11::shared_ptr;


D2K_NAMESPACE_OPEN


/** Demangle c++ names. */
std::string demangle(const char *name);

/**
 * This function collects some time utilities.
 *
 *  All measures are stored in seconds.
 *  Usage: get_start_time() should be used before get_end_time() and viceversa.
 */
class TimeUtilities
{
public:

  TimeUtilities()
    :
    status(true),
    times()
  {}

  /**
   * It freezes the thread for t milliseconds.
   */
  void sleep(unsigned int t);

  /**
   * It sets the start time for a measure.
   */
  void get_start_time();

  /**
   * It sets the end time for a measure.
   */
  void get_end_time();

  /**
   * It returns the number of measures done
   */
  int get_num_measures();

  /**
   * An overload of the operator [] is provided to access to all measures.
   */
  double &operator[] (const unsigned int num)
  {
    AssertThrow( num < times.size(),
                 ExcMessage("Invalind number num. It is higher than the number of times.") );

    return times[num];
  };

private:
  bool status; // This is used to check whether we are taking a start time or
  // and end time.
  std::chrono::high_resolution_clock::time_point t_start;
  std::chrono::high_resolution_clock::time_point t_end;
  std::vector<double> times;
};

/**
 * This function copy the text contained in @p in_file to the file
 * @p out_file .
 */
void append_to_file(const std::string &in_file, const std::string &out_file);

/**
 * A function that return the index of the first non existing folder matching
 * a pattern make by @p base and @p n_digits number. (base000, base001, base002, ...)
 * The research of the index starts from the value @p start and ends when @p index_max
 * is reached.
 */
unsigned int get_next_available_index_directory_name(const std::string &base, int n_digits=3,
                                                     unsigned int start=0, unsigned int index_max = 1000);

/**
 * A function that return the name of the first non existing folder matching
 * a pattern make by @p base and @p n_digits number. (base000, base001, base002, ...)
 * The research of the index starts from the value @p start and ends when @p index_max
 * is reached.
 */
std::string get_next_available_directory_name(const std::string &base, int n_digits=3,
                                              unsigned int start=0, unsigned int index_max = 1000);

/**
 * A function to check the existence of @p dir directory.
 */
bool dir_exists(const std::string &dir);

/**
 * A function to check the existence of @p file file.
 */
bool file_exists(const std::string &file);

/**
 * A function to create directory. It creates all directories needed.
 */
bool create_directory(const std::string &name);

/**
 * A function to copy a list of @p file ( "file1 file2 file3" ) in the
 * destination folder (@p destination)
 */
bool copy_files(const std::string &files, const std::string &destination);

/**
 * A function to make a copy of @p file with the name @p destination
 */
bool copy_file(const std::string &files, const std::string &destination);

/**
 * A function to rename a @p file with a new name @p new_file
 */
bool rename_file(const std::string &file, const std::string &new_file);

/// Cannot execute std::system(command)
DeclException1(ExcCannottExecuteCommand, std::string,
               << "Cannot execute command " << arg1 << "\n Please verify you have "
               << "the needed permissions.");

// Forward declaration for OverWriteStream:
template<typename Stream = std::ostream> class OverWriteStream;
/**
 * This class uses @p n_lines lines of @p stream_out to show the output.
 * Everytime it reaches the last line it comes back to the first line and
 * rewrites the line.
 *
 * The constructor takes as argument a stream for the output @p stream_out
 * (default = std::cout), the number of lines @p n_lines, and the width of the
 * output.
 *
 */
template<typename Stream>
class OverWriteStream
{
public:
  OverWriteStream( unsigned int n_lines = 1,
                   Stream &stream_out = std::cout,
                   unsigned int width = 60)
    :
    n_lines(n_lines),
    width(width),
    current_line(0),
    clear_next(true),
    stream_out(stream_out)
  {
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    rows_shell = w.ws_row;
    cols_shell = w.ws_col;
    stream_out.width(width);
  }

  /**
   * Flush output at the end so that we don't make a mess with the
   * console.
   */
  ~OverWriteStream()
  {
    for (; current_line<n_lines; ++current_line)
      stream_out << std::endl;
  }

  template<typename OBJ>
  OverWriteStream<Stream> &operator<<(OBJ &o)
  {
    clear();
    stream_out << o;
    return *this;
  }

  OverWriteStream<Stream> &operator<< (std::ostream& (*p) (std::ostream &))
  {
    class QueryStreambuf : public std::streambuf
    {
      // Implement a minimalistic stream buffer that only stores the fact
      // whether overflow or sync was called
    public:
      QueryStreambuf()
        : flushed_(false), newline_written_(false)
      {
      }
      bool flushed()
      {
        return flushed_;
      }
      bool newline_written()
      {
        return newline_written_;
      }
    private:
      int_type overflow(int_type ch)
      {
        newline_written_ = true;
        return ch;
      }
      int sync()
      {
        flushed_ = true;
        return 0;
      }

      bool flushed_;
      bool newline_written_;
    } query_streambuf;

    {
      // and initialize an ostream with this streambuf:
      std::ostream inject (&query_streambuf);
      inject << p;
    }

    if (query_streambuf.newline_written())
      current_line++;

    if (current_line == n_lines)
      {
        current_line = 0;
        for (unsigned int i=0; i<n_lines; ++i)
          stream_out << "\033[A"  << "\r";
        clear_next = true;
      }

    stream_out << p;

    return *this;
  }

  template <typename S, typename T> friend OverWriteStream<S>
  &operator << (OverWriteStream<S> &, const T &);

  Stream &get_stream()
  {
    return stream_out;
  }

  /**
   * Move the cursor at the end of the output. It is needed to avoid to rewrite
   * useful lines.
   * Finally, this method delete the class.
   */
  void end()
  {
    while (current_line <= n_lines-1)
      {
        stream_out << std::endl;
        current_line++;
      };
  }

  /**
   * Clear the next n lines, return back to the original point, and
   * reset the clear flag if force is set to true, otherwise do this
   * only if the internal counter is set to zero, i.e., we are at the
   * beginning of the next n lines.
   */
  void clear(bool force=false)
  {
    if (force)
      {
        while (current_line < n_lines-1)
          {
            stream_out << "\033[B";
            current_line++;
          };

        if (current_line>0)
          {
            for (unsigned int i=0; i<n_lines; ++i)
              {
                stream_out << "\r" << std::setfill(' ') << std::setw(width) << " " << "\r" << "\033[A";
                current_line--;
                if (current_line==0)
                  break;
              }
            stream_out << "\r" << std::setfill(' ') << std::setw(width) << " " << "\r";
          }
        else
          {
            stream_out << "\r" << std::setfill(' ') << std::setw(width) << " " << "\r";
          }
        clear_next = false;
      }
    else if (clear_next)
      {
        stream_out << "\r" << std::setfill(' ') << std::setw(width) << " " << "\r";
        stream_out << "\033[A" << std::endl;
        clear_next = false;
      }
  }

  /**
   * It returns the number of rows of the current shell.
   */
  unsigned int get_shell_rows()
  {
    return rows_shell;
  }

  /**
   * It returns the number of coloumns of the current shell.
   */
  unsigned int get_shell_cols()
  {
    return cols_shell;
  }

  /**
   * It returns current line of the stream.
   */
  int get_current_line()
  {
    return current_line;
  }

private:
  unsigned int cols_shell;
  unsigned int rows_shell;
  // total number of lines:
  const unsigned int n_lines;
  const unsigned int width;
  bool clear_next;
  // the current line:
  unsigned int current_line;
  // stream where the output will be written
  Stream &stream_out;

};

template <typename S, typename T>
inline
OverWriteStream<S> &operator<< (OverWriteStream<S> &output_stream, const T &t)
{
  output_stream.clear();
  output_stream.get_stream() << t;
  return output_stream;
}

// ======================================================================
// Explicit template functions. Only present in the include file.
// ======================================================================

/**
 * A simple function that accepts a vector as an input and returns a
 * second vector containing only the unique value among consecutive entries
 * of the original vector.
 */
template<class T>
std::vector<T> unique(const std::vector<T> &myvector)
{
  std::vector<T> ret;
  std::unique_copy(myvector.begin(), myvector.end(), std::back_inserter(ret));
  return ret;
}


/**
 * Return a string containing the content of the vector, with elements
 * separated by the @ sep parameter.
 */
template<class Type>
std::string print(const std::vector<Type> &list, const std::string sep=",")
{
  std::stringstream ret;
  if (list.size() > 0)
    ret << list[0];

  for (unsigned int i=1; i<list.size(); ++i)
    ret << sep << list[i];

  return ret.str();
}


/**
 * Return a string containing the content of the Point, with elements
 * separated by the @ sep parameter.
 */
template<int dim>
std::string print(const Point<dim> &point, const std::string sep=",")
{
  std::stringstream ret;
  ret << point[0];

  for (unsigned int i=1; i<dim; ++i)
    ret << sep << point[i];

  return ret.str();
}




/**
 * Return a human readable name of the type passed as argument.
 */
template <class T>
std::string type(const T &t)
{
  return demangle(typeid(t).name());
}

/**
 *  Construct a shared pointer to a non const class T. This is a
 *  convenience function to simplify the construction of shared
 *  pointers (which should replace dealii::SmartPointers):
 *
 *  @code
 *
 *  std_cxx11::shared_ptr<MyClass> my_ptr;
 *
 *  ...
 *
 *  my_ptr = SP(new MyClass);
 *
 *  @endcode
 */
template <class T>
inline shared_ptr<T>
SP(T *t)
{
  return shared_ptr<T>(t);
}

/**
 *  Construct a shared pointer to a const class T. This is a
 *  convenience function to simplify the construction of shared
 *  pointers (which should replace dealii::SmartPointers):
 *
 *  @code
 *
 *  std_cxx11::shared_ptr<const MyClass> my_ptr;
 *
 *  ...
 *  const MyClass * p = new MyClass;
 *  my_ptr = SP(p);
 *
 *  @endcode
 */
template <class T>
inline shared_ptr<const T>
SP(const T *t)
{
  return shared_ptr<const T>(t);
}

/**
 *  A simple class to shift a vector by a scalar.
 *  This function is deprecated in deal but needed in many codes
 */

template <typename VEC>
void vector_shift(VEC &in_vec, double a_scalar)
{
  for (auto i : in_vec.locally_owned_elements())
    in_vec[i] += a_scalar;
}

#ifdef D2K_WITH_SUNDIALS

#ifdef DEAL_II_WITH_MPI

#ifdef DEAL_II_WITH_TRILINOS
void copy(TrilinosWrappers::MPI::Vector &dst, const N_Vector &src);
void copy(N_Vector &dst, const TrilinosWrappers::MPI::Vector &src);
void copy(TrilinosWrappers::MPI::BlockVector &dst, const N_Vector &src);
void copy(N_Vector &dst, const TrilinosWrappers::MPI::BlockVector &src);
#endif // DEAL_II_WITH_TRILINOS

#ifdef DEAL_II_WITH_PETSC
void copy(PETScWrappers::MPI::Vector &dst, const N_Vector &src);
void copy(N_Vector &dst, const PETScWrappers::MPI::Vector &src);
void copy(PETScWrappers::MPI::BlockVector &dst, const N_Vector &src);
void copy(N_Vector &dst, const PETScWrappers::MPI::BlockVector &src);
#endif // DEAL_II_WITH_PETSC

#endif

void copy(BlockVector<double> &dst, const N_Vector &src);
void copy(N_Vector &dst, const BlockVector<double> &src);

#endif

#ifdef DEAL_II_WITH_TRILINOS
/**
 * A mixed deal.II - Trilinos monitor.. You can instantiate one object
 * of this type, and then call, in each function you want to monitor, the
 * method `auto t = timer.scoped_timer("Section");` which will automatically
 * start the timer "Section", and stops it when the object `t` is destroyed,
 * using both a dealii::TimerOutput and a Teuchos::TimeMonitor.
 */
class TimeMonitor
{
public:
  /**
   * A mixed deal.II - Trilinos monitor.
   */
  TimeMonitor(const MPI_Comm &comm=MPI_COMM_WORLD,
              std::ostream  &stream=std::cout) :
    outstream(stream),
    out(outstream, Utilities::MPI::this_mpi_process(comm) == 0),
    dealii_timer(comm, out, TimerOutput::never, TimerOutput::cpu_and_wall_times)
  {}

  /**
   * Print both deal.II and Trilinos Summaries.
   */
  ~TimeMonitor()
  {
    dealii_timer.print_summary();
    Teuchos::TimeMonitor::summarize(outstream);
  }

  /**
   * Helper class to enter/exit sections. The only role of this class is to
   * create an object that upon construction will start the given timers,
   * and upon destruction will stop them.
   */
  class Scope
  {
  public:
    /**
     * Start Trilinos and deal.II scoped monitors.
     */
    Scope(Teuchos::Time &trilinos_time,
          TimerOutput &dealii_timer_output,
          const std::string &section) :
      trilinos_scoped_monitor(trilinos_time),
      dealii_scoped_monitor(dealii_timer_output, section)
    {}

  private:
    /**
     * Trilinos version of scoped monitor.
     */
    Teuchos::TimeMonitor trilinos_scoped_monitor;

    /**
     * Deal.II version of scoped monitor.
     */
    TimerOutput::Scope dealii_scoped_monitor;
  };

  /**
   * Create and start a timer named after the parameter @p section. When the
   * created object is destroyed, the timer is stopped. Repeated calls with
   * the same timer are accumulated, together with some statistics.
   */
  Scope
  scoped_timer(const std::string &section) const
  {
    if (timers.find(section) == timers.end())
      timers[section] = Teuchos::TimeMonitor::getNewCounter(section);

    return Scope(*(timers[section]), dealii_timer, section);
  }

private:
  /**
   * Output stream.
   */
  std::ostream &outstream;

  /**
   * Conditional output stream. Output is only generated on processor 0.
   */
  ConditionalOStream out;

  /**
   * A list of Trilinos timers.
   */
  mutable std::map<std::string, Teuchos::RCP<Teuchos::Time> > timers;

  /**
   * The deal.II TimeMonitor object.
   */
  mutable TimerOutput dealii_timer;
};
#endif

D2K_NAMESPACE_CLOSE

#endif

