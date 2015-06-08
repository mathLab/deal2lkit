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
 * A function to check the existence of @p file
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

#ifdef DEAL_II_WITH_TRILINOS

/**
 *  Define namespace "SacadoUtilities"
 *  some utilities aimed at solving non-linear problems
 *  are implemented.
 *
 */

#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> Sdouble;
typedef Sacado::Fad::DFad<Sdouble> SSdouble;

namespace SacadoUtilities
{
  /**
   *  Extract local dofs values and initialize the number
   *  of independent variables up to the second order derivative.
   *
   */
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
        else if (typeid(Number) == typeid(Sdouble))
          {
            Sdouble ildv = global_vector (local_dof_indices[i]);
            ildv.diff (i, dofs_per_cell);
            ((Sdouble &)independent_local_dofs[i]) = ildv;
          }
        else if (typeid(Number) == typeid(SSdouble))
          {
            Sdouble ildv = global_vector (local_dof_indices[i]);
            ildv.diff (i, dofs_per_cell);
            ((SSdouble &)independent_local_dofs[i]) = ildv;
            ((SSdouble &)independent_local_dofs[i]).diff(i,dofs_per_cell);
          }
        else
          {
            Assert (false, ExcNotImplemented());
          }
      }
  }

  /**
   *  Extract independent local dofs values of a vector_variable
   *  in a cell and store them in
   *  std::vector <Tensor <1, spacedim, Number> > us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_values (const FEValues<dim, spacedim> &fe_values,
              const std::vector<Number> &independent_local_dof_values,
              std::vector <Tensor <1, spacedim, Number> > &us,
              const FEValuesExtractors::Vector &vector_variable)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            us[q] += independent_local_dof_values[i]*fe_values[vector_variable].value(i,q);
          }
      }
  }

  /**
   *  Extract independent local dofs values of a scalar_variable
   *  in a cell and store them in
   *  std::vector <Number>  us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_values (const FEValues<dim, spacedim> &fe_values,
              const std::vector<Number> &independent_local_dof_values,
              std::vector <Number> &us,
              const FEValuesExtractors::Scalar &scalar_variable)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            us[q] += independent_local_dof_values[i]*fe_values[scalar_variable].value(i,q);
          }
      }
  }

  /**
   *  Extract independent local dofs values of a vector_variable
   *  on face of a cell and store them in
   *  std::vector <Tensor <1, spacedim, Number> > us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_face_values (const FEFaceValues<dim, spacedim> &fe_values,
                   const std::vector<Number> &independent_local_dof_values,
                   std::vector <Tensor <1, spacedim, Number> > &us,
                   const FEValuesExtractors::Vector &vector_variable)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            us[q] += independent_local_dof_values[i]*fe_values[vector_variable].value(i,q);
          }
      }
  }

  /**
   *  Extract independent local dofs values of a scalar_variable
   *  on face of a cell and store them in
   *  std::vector <Number>  us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_face_values (const FEFaceValues<dim, spacedim> &fe_values,
                   const std::vector<Number> &independent_local_dof_values,
                   std::vector <Number> &us,
                   const FEValuesExtractors::Scalar &scalar_variable)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            us[q] += independent_local_dof_values[i]*fe_values[scalar_variable].value(i,q);
          }
      }
  }


  /**
   *  Extract divergence values of a vector_variable
   *  on face of a cell and store them in
   *  std::vector <Number>  us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_div_values (const FEValues<dim, spacedim> &fe_values,
                  const std::vector<Number> &independent_local_dof_values,
                  std::vector <Number>  &us,
                  const FEValuesExtractors::Vector &vector_variable)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            us[q] += independent_local_dof_values[i]*fe_values[vector_variable].divergence(i,q);
          }
      }
  }


  /**
   *  Extract gradient values of a vector_variable
   *  on face of a cell and store them in
   *  std::vector <Tensor <2, spacedim, Number> > grad_us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_grad_values (const FEValues<dim, spacedim> &fe_values,
                   const std::vector<Number> &independent_local_dof_values,
                   std::vector <Tensor <2, spacedim, Number> > &grad_us,
                   const FEValuesExtractors::Vector &vector_variable)
  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(grad_us.size(), n_q_points);

    std::vector<Tensor<2,spacedim,double> > buffer(n_q_points);

    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        buffer[q] += fe_values[vector_variable].gradient(i,q);

    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned d=0; d<spacedim; ++d)
        for (unsigned dd=0; dd<spacedim; ++dd)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            grad_us[q][d][dd] += independent_local_dof_values[i]*buffer[q][d][dd];
  }

  /**
   *  Compute the deformation gradient F in each quadrature point
   *  of a cell and it is stored in
   *  std::vector <Tensor <2, spacedim, Number> > Fs
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_F_values (const FEValues<dim, spacedim> &fe_values,
                const std::vector<Number> &independent_local_dof_values,
                std::vector <Tensor <2, spacedim, Number> > &Fs,
                const FEValuesExtractors::Vector &vector_variable)
  {
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(Fs.size(), n_q_points);

    SacadoUtilities::get_grad_values (fe_values, independent_local_dof_values, Fs, vector_variable);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int d=0; d<dim; ++d)
          Fs[q][d][d] += 1.0; // I + grad(u)
      }
  }

  /**
   *  Extract symmetric gradient values of a vector_variable
   *  on face of a cell and store them in
   *  std::vector <Tensor <2, spacedim, Number> > grad_us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_sym_grad_values (const FEValues<dim, spacedim> &fe_values,
                       const std::vector<Number> &independent_local_dof_values,
                       std::vector <Tensor <2, spacedim, Number> > &grad_us,
                       const FEValuesExtractors::Vector &vector_variable)
  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(grad_us.size(), n_q_points);

    std::vector<Tensor<2,spacedim,double> > buffer(n_q_points);

    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        buffer[q] += fe_values[vector_variable].symmetric_gradient(i,q);

    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned d=0; d<spacedim; ++d)
        for (unsigned dd=0; dd<spacedim; ++dd)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            grad_us[q][d][dd] += independent_local_dof_values[i]*buffer[q][d][dd];
  }
}// end namespace
#endif // DEAL_II_WITH_TRILINOS


#endif
