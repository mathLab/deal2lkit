#ifndef UTILITIES_HLT
#define UTILITIES_HLT

#include <deal.II/base/utilities.h>
#include <deal.II/base/smartpointer.h>
#include <typeinfo>
#include <cxxabi.h>

using namespace dealii;


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

// Anonymous namespace, to hide implementation detail for the type
// function below
namespace
{
  struct handle
  {
    char *p;
    handle(char *ptr) : p(ptr) { }
    ~handle()
    {
      delete p;
    }
  };

  std::string demangle(const char *name)
  {
    int status = -4; // some arbitrary value to eliminate the compiler warning
    handle result( abi::__cxa_demangle(name, NULL, NULL, &status) );
    return (status==0) ? result.p : name ;
  }
}

/**
 * Return a human readable name of the type passed as argument.
 */
template <class T>
std::string type(const T &t)
{
  return demangle(typeid(t).name());
}

/** Construct a unique pointer to a non const class T. */
template <class T>
std::unique_ptr<T>
UP(T *t)
{
  return std::unique_ptr<T>(t);
}

/** Construct a unique pointer to a const class T. */
template <class T>
std::unique_ptr<const T>
UP(const T *t)
{
  return std::unique_ptr<const T>(t);
}

/** Construct a shared pointer to a non const class T. */
template <class T>
std::shared_ptr<T>
SP(T *t)
{
  return std::shared_ptr<T>(t);
}

/** Construct a shared pointer to a const class T. */
template <class T>
std::shared_ptr<const T>
SP(const T *t)
{
  return std::shared_ptr<const T>(t);
}


#endif
