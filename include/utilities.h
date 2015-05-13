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

#ifdef DEAL_II_WITH_TRILINOS

/**
 *
 *  Extract local dofs values and initialize the number
 *  of independent variables up to the second order derivative.
 *
 *  @begin code
 *
 *  ...
 *
 *  std::vector<types::global_dof_index>    local_dof_indices (dofs_per_cell);
 *  std::vector<Number> independent_local_dof_values (dofs_per_cell);
 *
 *  fe_values.reinit (cell);
 *  cell->get_dof_indices (local_dof_indices);
 *
 *  extract_local_dofs (solution, local_dof_indices, independent_local_dof_values);
 *
 *  ...
 *
 *  @end
 */

#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> Sdouble;
typedef Sacado::Fad::DFad<Sdouble> SSdouble;

template <typename Number, typename VEC>
void
extract_local_dofs (const VEC &global_vector,
                    const std::vector<types::global_dof_index> &local_dof_indices,
                    std::vector<Number> &independent_local_dofs)
{
  const unsigned int dofs_per_cell = local_dof_indices.size();
  for (unsigned int i=0; i < dofs_per_cell; ++i)
    {
      if (typeid(Number) == typeid(double))
        {
          independent_local_dofs[i] = global_vector (local_dof_indices[i]);
        }
      else
        {
          Sdouble ildv = global_vector (local_dof_indices[i]);
          ildv.diff (i, dofs_per_cell);
          ((Sdouble &)independent_local_dofs[i]) = ildv;
          if (typeid(Number) == typeid(SSdouble))
            {
              ((SSdouble &)independent_local_dofs[i]).diff(i,dofs_per_cell);
            }
        }
    }
}

#endif // DEAL_II_WITH_TRILINOS


#endif
