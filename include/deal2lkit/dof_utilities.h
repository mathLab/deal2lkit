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

#ifndef _d2k_dof_utilities_h
#define _d2k_dof_utilities_h

#include <deal2lkit/config.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/smartpointer.h>
#include <typeinfo>
#include <deal.II/base/exceptions.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>


#include <sstream>

using namespace dealii;



/**
 *  Define namespace "DOFUtilities"
 *  some utilities aimed at solving non-linear problems
 *  are implemented.
 *
 */

#ifdef DEAL_II_WITH_TRILINOS


#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> Sdouble;
typedef Sacado::Fad::DFad<Sdouble> SSdouble;
#endif

D2K_NAMESPACE_OPEN
namespace DOFUtilities
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
                      std::vector<Number> &independent_local_dofs,
                      const unsigned int n_vectors=1,
                      const unsigned int my_index=0)
  {
    AssertDimension(local_dof_indices.size(),independent_local_dofs.size());

    const unsigned int dofs_per_cell = local_dof_indices.size();

    for (unsigned int i=0; i < dofs_per_cell; ++i)
      {
        if (typeid(Number) == typeid(double))
          {
            independent_local_dofs[i] = global_vector (local_dof_indices[i]);
          }
#ifdef DEAL_II_WITH_TRILINOS
        else if (typeid(Number) == typeid(Sdouble))
          {
            Sdouble ildv(n_vectors*dofs_per_cell, i+my_index*dofs_per_cell,
                         global_vector (local_dof_indices[i]));
            ((Sdouble &)independent_local_dofs[i]) = ildv;
          }
        else if (typeid(Number) == typeid(SSdouble))
          {
            SSdouble ildv(n_vectors*dofs_per_cell, i+my_index*dofs_per_cell,
                          global_vector(local_dof_indices[i]));
            ildv.val() = Sdouble(n_vectors*dofs_per_cell, i+my_index*dofs_per_cell,
                                 global_vector(local_dof_indices[i]));
            ((SSdouble &)independent_local_dofs[i]) = ildv;
          }
#endif
        else
          {
            Assert (false, ExcNotImplemented());
          }
      }
  }

  /**
   *  Extract independent local dofs values
   *  in a cell and store them in
   *  std::vector <std::vector<Number> > us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_values (const FEValuesBase<dim, spacedim> &fe_values,
              const std::vector<Number> &independent_local_dof_values,
              std::vector <std::vector<Number> > &us)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;
    const unsigned int           n_components  = fe_values.get_fe().n_components();

    AssertDimension(us.size(), n_q_points);
    AssertDimension(us[0].size(), n_components);
    AssertDimension(independent_local_dof_values.size(), dofs_per_cell);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        us[q] = std::vector<Number> (n_components,0.0);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<n_components; ++j)
            {
              FEValuesExtractors::Scalar s(j);
              us[q][j] += independent_local_dof_values[i]*fe_values[s].value(i,q);
            }
      }
  }

  /**
   *  Extract independent local dofs values of a vector_variable
   *  in a cell and store them in
   *  std::vector<Tensor<1,spacedim,Number>> us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>

  void
  get_values (const FEValuesBase<dim, spacedim> &fe_values,
              const std::vector<Number> &independent_local_dof_values,
              const FEValuesExtractors::Vector &vector_variable,
              std::vector <Tensor <1, spacedim, Number> > &us)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);
    AssertDimension(independent_local_dof_values.size(), dofs_per_cell);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        us[q] = 0;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          if (fe_values[vector_variable].value(i,q).norm() > 0.0)
            for (unsigned int j=0; j<spacedim; ++j)
              us[q][j] += independent_local_dof_values[i]*fe_values[vector_variable].value(i,q)[j];
      }
  }

  /**
   *  Extract independent local dofs values of a scalar_variable
   *  in a cell and store them in
   *  std::vector<Number>  us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_values (const FEValuesBase<dim, spacedim> &fe_values,
              const std::vector<Number> &independent_local_dof_values,
              const FEValuesExtractors::Scalar &scalar_variable,
              std::vector <Number> &us)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);
    AssertDimension(independent_local_dof_values.size(), dofs_per_cell);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        us[q] = 0;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          us[q] += independent_local_dof_values[i]*fe_values[scalar_variable].value(i,q);
      }
  }


  /**
   *  Extract divergence values of a vector_variable
   *  on face of a cell and store them in
   *  std::vector<Number> us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_divergences (const FEValuesBase<dim, spacedim> &fe_values,
                   const std::vector<Number> &independent_local_dof_values,
                   const FEValuesExtractors::Vector &vector_variable,
                   std::vector <Number>  &us)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);
    AssertDimension(independent_local_dof_values.size(), dofs_per_cell);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        us[q] = 0;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          us[q] += independent_local_dof_values[i]*fe_values[vector_variable].divergence(i,q);
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
  get_gradients (const FEValuesBase<dim, spacedim> &fe_values,
                 const std::vector<Number> &independent_local_dof_values,
                 const FEValuesExtractors::Vector &vector_variable,
                 std::vector <Tensor <2, spacedim, Number> > &grad_us)
  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(grad_us.size(), n_q_points);
    AssertDimension(independent_local_dof_values.size(), dofs_per_cell);


    for (unsigned int q=0; q<n_q_points; ++q)
      {
        grad_us[q] = 0;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int d=0; d<spacedim; ++d)
            for (unsigned int dd=0; dd<spacedim; ++dd)
              grad_us[q][d][dd] += independent_local_dof_values[i]*fe_values[vector_variable].gradient(i,q)[d][dd];
      }

  }


  /**
   *  Extract curl values of a vector_variable
   *  on face or cell and store them in
   *  std::vector <Tensor <1, spacedim, Number> > curl_us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_curls (const FEValuesBase<dim, spacedim> &fe_values,
             const std::vector<Number> &independent_local_dof_values,
             const FEValuesExtractors::Vector &vector_variable,
             std::vector <Tensor <1, (spacedim > 2 ? spacedim : 1), Number> > &curl_us)
  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(curl_us.size(), n_q_points);
    AssertDimension(independent_local_dof_values.size(), dofs_per_cell);


    for (unsigned int q=0; q<n_q_points; ++q)
      {
        curl_us[q] = 0;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            auto c = fe_values[vector_variable].curl(i,q);
            for (unsigned int d=0; d<(spacedim > 2 ? spacedim : 1); ++d)
              curl_us[q][d] += independent_local_dof_values[i]*c[d];
          }
      }

  }


  /**
   *  Extract gradient values of a scalar_variable
   *  on face of a cell and store them in
   *  std::vector <Tensor <1, spacedim, Number> > grad_us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_gradients (const FEValuesBase<dim, spacedim> &fe_values,
                 const std::vector<Number> &independent_local_dof_values,
                 const FEValuesExtractors::Scalar &scalar_variable,
                 std::vector <Tensor <1, spacedim, Number> > &grad_us)
  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(grad_us.size(), n_q_points);
    AssertDimension(independent_local_dof_values.size(), dofs_per_cell);


    for (unsigned int q=0; q<n_q_points; ++q)
      {
        grad_us[q] = 0;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int d=0; d<spacedim; ++d)
            grad_us[q][d] += independent_local_dof_values[i]*fe_values[scalar_variable].gradient(i,q)[d];
      }

  }


  /**
   *  Compute the deformation gradient in each quadrature point
   *  of a cell and it is stored in
   *  std::vector <Tensor <2, spacedim, Number> > Fs
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>
  void
  get_deformation_gradients (const FEValuesBase<dim, spacedim> &fe_values,
                             const std::vector<Number> &independent_local_dof_values,
                             const FEValuesExtractors::Vector &vector_variable,
                             std::vector <Tensor <2, spacedim, Number> > &Fs)
  {
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(Fs.size(), n_q_points);

    DOFUtilities::get_gradients (fe_values, independent_local_dof_values, vector_variable, Fs);

    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int d=0; d<dim; ++d)
        Fs[q][d][d] += 1.0; // I + grad(u)
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
  get_symmetric_gradients (const FEValuesBase<dim, spacedim> &fe_values,
                           const std::vector<Number> &independent_local_dof_values,
                           const FEValuesExtractors::Vector &vector_variable,
                           std::vector <Tensor <2, spacedim, Number> > &grad_us)
  {
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(grad_us.size(), n_q_points);

    DOFUtilities::get_gradients (fe_values, independent_local_dof_values, vector_variable, grad_us);
    for (unsigned int q=0; q<n_q_points; ++q)
      {
        grad_us[q] += transpose(grad_us[q]);
        grad_us[q] *= 0.5;
      }

  }


  /**
   *  Extract laplacians local dofs values of a scalar variable
   *  in a cell and store them in
   *  std::vector<Number> us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>

  void
  get_laplacians (const FEValuesBase<dim, spacedim> &fe_values,
                  const std::vector<Number> &independent_local_dof_values,
                  const FEValuesExtractors::Scalar &variable,
                  std::vector <Number> &us)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);
    AssertDimension(independent_local_dof_values.size(), dofs_per_cell);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        us[q] = 0;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          us[q] += independent_local_dof_values[i]*trace(fe_values[variable].hessian(i,q));
      }
  }


  /**
   *  Extract laplacians local dofs values of a vector variable
   *  in a cell and store them in
   *  std::vector <Tensor<1,Number > > us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>

  void
  get_laplacians (const FEValuesBase<dim, spacedim> &fe_values,
                  const std::vector<Number> &independent_local_dof_values,
                  const FEValuesExtractors::Vector &variable,
                  std::vector <Tensor<1,spacedim,Number> > &us)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);
    AssertDimension(independent_local_dof_values.size(), dofs_per_cell);

    for (unsigned int d=0; d<spacedim; ++d)
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          us[q][d] = 0;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            us[q][d] += independent_local_dof_values[i]*trace(fe_values[variable].hessian(i,q)[d]);
        }
  }

  /**
   *  Extract hessians local dofs values of a scalar variable
   *  in a cell and store them in
   *  std::vector <Tensor<2,Number > > us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>

  void
  get_hessians (const FEValuesBase<dim, spacedim> &fe_values,
                const std::vector<Number> &independent_local_dof_values,
                const FEValuesExtractors::Scalar &variable,
                std::vector <Tensor<2,spacedim,Number> > &us)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);
    AssertDimension(independent_local_dof_values.size(), dofs_per_cell);

    for (unsigned int d=0; d<spacedim; ++d)
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          us[q]= 0;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            us[q] += independent_local_dof_values[i]*fe_values[variable].hessian(i,q);
        }
  }


  /**
   *  Extract hessians local dofs values of a vector variable
   *  in a cell and store them in
   *  std::vector <Tensor<3,Number > > us
   *  whose size is number of quadrature points in the cell.
   *
   */
  template <int dim, int spacedim, typename Number>

  void
  get_hessians (const FEValuesBase<dim, spacedim> &fe_values,
                const std::vector<Number> &independent_local_dof_values,
                const FEValuesExtractors::Vector &variable,
                std::vector <Tensor<3,spacedim,Number> > &us)

  {
    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    AssertDimension(us.size(), n_q_points);
    AssertDimension(independent_local_dof_values.size(), dofs_per_cell);

    for (unsigned int d=0; d<spacedim; ++d)
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          us[q]= 0;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            us[q] += independent_local_dof_values[i]*fe_values[variable].hessian(i,q);
        }
  }



}// end namespace


D2K_NAMESPACE_CLOSE

#endif
