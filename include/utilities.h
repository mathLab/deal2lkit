#ifndef UTILITIES_HLT
#define UTILITIES_HLT

#include <deal.II/base/utilities.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/fe/fe.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>

#include <typeinfo>
#include <cxxabi.h>
#include <sstream>
#include <sys/ioctl.h>

using namespace dealii;
using std_cxx11::shared_ptr;

/**
 * SmartPointers are usually not used to point to objects created with
 * new. However, sometimes this is useful. The distruction of a
 * SmartPointer requires to split the step in two parts. This little
 * utility does precisely this.
 *
 * @deprecated SmartPointers have been supersed by
 * std_cxx11::shared_ptr, which takes care of destruction as well.
 */
template <typename TYPE>
void smart_delete (SmartPointer<TYPE> &sp) DEAL_II_DEPRECATED;

/** Demangle c++ names. */
std::string demangle(const char *name);

/**
 * This function collects some time utilities:
 *  - sleep(unsigned int t) freezes the thread for t milliseconds;
 *  - get_start_time() sets the start time for a measure
 *  - get_end_time() sets the end time for a measure
 *  - get_num_measures() return the number of measures done
 *  - an overload of the operator [] is provided to access to all measures
 *
 *  All measures are stored in seconds.
 *  get_start_time() should be used before get_end_time() and viceversa.
 */
class TimeUtilities
{
public:

  TimeUtilities()
    :
    status(true),
    times()
  {}

  void sleep(unsigned int t);
  void get_start_time();
  void get_end_time();
  int get_num_measures();
  double &operator[] (const int num)
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
 * This function copyt the text contained in @p in_file to the file
 * @p out_file .
 */
void append_to_file(const std::string &in_file, const std::string &out_file);

/**
 * A function that return the index of the first non existing folder matching
 * a pattern make by @p base and @p n_digits number. (base000, base001, base002, ...)
 */
int get_next_available_index_directory_name(const std::string &base, int n_digits=3);

/**
 * A function that return the name of the first non existing folder matching
 * a pattern make by @p base and @p n_digits number. (base000, base001, base002, ...)
 */
std::string get_next_available_directory_name(const std::string &base, int n_digits=3);

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

/**
 * This class rewrite @p n_max lines of output
 */
template<typename Stream>
class FilteredStream
{
public:
  FilteredStream( Stream &stream_out = std::cout,
                  unsigned int n_lines = 1,
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
  };

  /**
   * Flush output at the end so that we don't make a mess with the
   * console.
   */
  ~FilteredStream()
  {
    for (; current_line<n_lines; ++current_line)
      stream_out << std::endl;
  };

  template<typename OBJ>
  FilteredStream<Stream> &operator<<(OBJ &o)
  {
    clear();
    stream_out << o;
    return *this;
  };

  FilteredStream<Stream> &operator<< (std::ostream& (*p) (std::ostream &))
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
          stream_out << "\e[A"  << "\r";
        clear_next = true;
      }

    stream_out << p;

    return *this;
  };

  template <typename S, typename T> friend FilteredStream<S>
  &operator << (FilteredStream<S> &, const T &);

  Stream &get_stream()
  {
    return stream_out;
  };

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
            stream_out << "\e[B";
            current_line++;
          }

        if (current_line>0)
          {
            for (unsigned int i=0; i<n_lines; ++i)
              {
                stream_out << "\r" << std::setfill(' ') << std::setw(width) << " " << "\r" << "\e[A";
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
        stream_out << "\e[A" << std::endl;
        clear_next = false;
      }
  };

  unsigned int get_shell_rows()
  {
    return rows_shell;
  }
  unsigned int get_shell_cols()
  {
    return cols_shell;
  }
  int get_current_line()
  {
    return current_line;
  }

private:
  unsigned int cols_shell;
  unsigned int rows_shell;
  // total number of lines:
  const unsigned int n_lines;
  // the current line:
  const unsigned int width;
  bool clear_next;
  unsigned int current_line;
  // stream where the output will be written
  Stream &stream_out;

};

template <typename S, typename T>
inline
FilteredStream<S> &operator<< (FilteredStream<S> &output_stream, const T &t)
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
 *  @begin code
 *
 *  std_cxx11::shared_ptr<MyClass> my_ptr;
 *
 *  ...
 *
 *  my_ptr = SP(new MyClass);
 *
 *  @end
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
 *  @begin code
 *
 *  std_cxx11::shared_ptr<const MyClass> my_ptr;
 *
 *  ...
 *  const MyClass * p = new MyClass;
 *  my_ptr = SP(p);
 *
 *  @end
 */
template <class T>
inline shared_ptr<const T>
SP(const T *t)
{
  return shared_ptr<const T>(t);
}

template <typename TYPE>
void smart_delete (SmartPointer<TYPE> &sp)
{
  if (sp)
    {
      TYPE *p = sp;
      sp = 0;
      delete p;
    }
}

#endif
