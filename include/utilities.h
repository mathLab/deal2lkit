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

int get_next_available_index_directory_name(const std::string &base, int n_digits=3);
std::string get_next_available_directory_name(const std::string &base, int n_digits=3);

bool exists_file(const std::string &file);

bool create_directory(const std::string &name);

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
   *  Extract independent local dofs values of a cell and store them in
   *  std::vector <Tensor <1, spacedim, Number> > us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_u_values (const FEValues<dim, spacedim> &fe_values,
                const std::vector<Number> &independent_local_dof_values,
                std::vector <Tensor <1, spacedim, Number> > &us)
  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int c = fe_values.get_fe().system_to_component_index(i).first;
            us[q][c] += independent_local_dof_values[i]*fe_values.shape_value_component(i,q,c);
          }
      }
  }

  /**
   *  Extract independent local dofs values of a cell,
   *  compute the gradient grad_us in each quadrature point,
   *  which is stored in
   *  std::vector <Tensor <2, spacedim, Number> > grad_us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_grad_u_values (const FEValues<dim, spacedim> &fe_values,
                     const std::vector<Number> &independent_local_dof_values,
                     std::vector <Tensor <2, spacedim, Number> > &grad_us)
  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(grad_us.size(), n_q_points);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int c = fe_values.get_fe().system_to_component_index(i).first;
            for (unsigned int d=0; d<dim; ++d)
              {
                grad_us[q][c][d] += independent_local_dof_values[i]*fe_values.shape_grad_component(i,q,c)[d]; // grad(u)
              }
          }
      }
  }

  /**
   *  Extract independent local dofs values of a cell,
   *  compute the deformation gradient F in each quadrature point,
   *  which is stored in
   *  std::vector <Tensor <2, spacedim, Number> > Fs
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_F_values (const FEValues<dim, spacedim> &fe_values,
                const std::vector<Number> &independent_local_dof_values,
                std::vector <Tensor <2, spacedim, Number> > &Fs)
  {
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(Fs.size(), n_q_points);

    SacadoUtilities::get_grad_u_values (fe_values, independent_local_dof_values, Fs);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int d=0; d<dim; ++d)
          Fs[q][d][d] += 1.0; // I + grad(u)
      }
  }

  /**
   *  Extract independent local dofs values of a face of a cell and store them in
   *  std::vector <Tensor <1, spacedim, Number> > us_face
   *  whose size is number of quadrature points in the cell face.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_u_face_values (const FEFaceValues<dim, spacedim> &fe_face_values,
                     const std::vector<Number> &independent_local_dof_values,
                     std::vector <Tensor <1, spacedim, Number> > &us_face)
  {
    const unsigned int           dofs_per_face = fe_face_values.dofs_per_cell;
    const unsigned int           n_face_q_points    = fe_face_values.n_quadrature_points;

    AssertDimension(us_face.size(), n_face_q_points);

    for (unsigned int qf=0; qf<n_face_q_points; ++qf)
      {
        for (unsigned int i=0; i<dofs_per_face; ++i)
          {
            const unsigned int c = fe_face_values.get_fe().system_to_component_index(i).first;
            us_face[qf][c] += independent_local_dof_values[i]*fe_face_values.shape_value_component(i,qf,c);
          }
      }
  }
}
#endif // DEAL_II_WITH_TRILINOS


#endif
