#ifndef UTILITIES_HLT
#define UTILITIES_HLT

#include <deal.II/base/utilities.h>
#include <deal.II/base/smartpointer.h>
#include <typeinfo>
#include <cxxabi.h>

#include <deal.II/base/std_cxx11/shared_ptr.h>

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

// ======================================================================
// Explicit template functions. Only present in the include file.
// ======================================================================

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
