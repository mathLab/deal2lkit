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
 */
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

/** Demangle c++ names. */
std::string demangle(const char *name);

/**
 * Return a human readable name of the type passed as argument.
 */
template <class T>
std::string type(const T &t)
{
  return demangle(typeid(t).name());
}

/** Construct a shared pointer to a non const class T. */
template <class T>
shared_ptr<T>
SP(T *t)
{
  return shared_ptr<T>(t);
}

/** Construct a shared pointer to a const class T. */
template <class T>
shared_ptr<const T>
SP(const T *t)
{
  return shared_ptr<const T>(t);
}


#endif
