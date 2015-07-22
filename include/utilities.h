#ifndef UTILITIES_HLT
#define UTILITIES_HLT

#include <deal.II/base/utilities.h>
#include <deal.II/base/smartpointer.h>
#include <typeinfo>
#include <cxxabi.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/fe/fe.h>

#include <deal.II/base/std_cxx11/shared_ptr.h>

#include <sstream>

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
