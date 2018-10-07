//-----------------------------------------------------------
//
//    Copyright (C) 2016 by the deal2lkit authors
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

#ifndef d2k_sacado_tools_h
#define d2k_sacado_tools_h

#include <deal2lkit/config.h>
#include <deal.II/base/config.h>


/**
 *  Define namespace "SacadoTools" some functions for working with
 *  Sacado-type variables.
 *
 */


#ifdef DEAL_II_WITH_TRILINOS


#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> Sdouble;
typedef Sacado::Fad::DFad<Sdouble> SSdouble;

#include <deal.II/base/tensor.h>

using namespace dealii;



D2K_NAMESPACE_OPEN

namespace SacadoTools
{
  /**
   * This function compute the scalar product
   * between two tensors with rank 1 of different types
   *  i.e., Number and double  returning a Number.
   */
  template<int dim, typename Number>
  Number scalar_product (const Tensor<1,dim,Number> &T1,
                         const Tensor<1,dim,double> &T2)
  {
    Number ret=0;
    for (unsigned int i=0; i<dim; ++i)
      ret += T1[i]*T2[i];
    return ret;
  }



  /**
   * This function compute the scalar product
   * between two tensors with rank 1 of different types
   *  i.e., Number and double  returning a Number.
   */
  template<int dim, typename Number>
  Number scalar_product (const Tensor<1,dim,double> &T1,
                         const Tensor<1,dim,Number> &T2)
  {
    return SacadoTools::scalar_product(T2,T1);
  }



  /**
   * This function is needed to avoid conflict
   */
  template<int dim>
  double scalar_product (const Tensor<1,dim,double> &T1,
                         const Tensor<1,dim,double> &T2)
  {
    return T1*T2;
  }



  /**
   * This function compute the scalar product
   * between two tensors with rank 2 of different types
   *  i.e., Number and double  returning a Number.
   */
  template<int dim, typename Number>
  Number scalar_product (const Tensor<2,dim,Number> &T1,
                         const Tensor<2,dim,double> &T2)
  {
    Number ret=0;
    for (unsigned int i=0; i<dim; ++i)
      for (unsigned int j=0; j<dim; ++j)
        ret += T1[i][j]*T2[i][j];
    return ret;
  }



  /**
   * This function compute the scalar product
   * between two tensors with rank 2 of different types
   *  i.e., Number and double  returning a Number.
   */
  template<int dim, typename Number>
  Number scalar_product (const Tensor<2,dim,double> &T1,
                         const Tensor<2,dim,Number> &T2)
  {
    return SacadoTools::scalar_product(T2,T1);
  }



  /**
   * This function is needed to avoid conflict
   */
  template<int dim>
  double scalar_product (const Tensor<2,dim,double> &T1,
                         const Tensor<2,dim,double> &T2)
  {
    return dealii::scalar_product(T1,T2);
  }



// Helper functions val() which recursively calls val() function of Sacado

  /**
   * Do nothing, required for compatibility
   */
  double val(const double &n)
  {
    return n;
  }



  /**
   * Call val on a Sacado<number> and return a number
   */
  template <typename number>
  number val(const Sacado::Fad::DFad<number> &n)
  {
    return n.val();
  }



  /**
    * Downgrade from Tensor<index,dim,Sacado<number> > to Tensor<index,dim,number>
    */
  template <int index, int dim, typename number>
  Tensor<index,dim,number> val(const Tensor<index,dim,Sacado::Fad::DFad<number> > &T)
  {
    Tensor<index,dim,number> ret;
    for (unsigned int i=0; i<dim; ++i)
      ret[i] = val(T[i]);
    return ret;
  }



  /**
   * Downgrade a std::vector<Tensor<index,dim,Sacado<number> > > to
   * std::vector<Tensor<index,dim, number> >
   */
  template <int index, int dim, typename number>
  std::vector<Tensor<index,dim,number> > val(const std::vector<Tensor<index,dim,Sacado::Fad::DFad<number> > > &T)
  {
    std::vector<Tensor<index,dim,number> > ret(T.size());
    for (unsigned int q=0; q<T.size(); ++q)
      ret[q] = val(T[q]);
    return ret;
  }



  /**
   * Downgrade a std::vector<Sacado<number> > to std::vector<number>
   */
  template <typename number>
  std::vector<number> val(const std::vector<Sacado::Fad::DFad<number> > v)
  {
    std::vector<number> ret(v.size());
    for (unsigned int i=0; i<v.size(); ++i)
      ret[i] = v[i].val();
    return ret;
  }



  // helper functions to_double wich downgrade Sacado type to double

  /**
   * Do nothing, required for compatibility.
   */
  double to_double(const double &s)
  {
    return s;
  }



  /**
   * Return a double given a Sdouble
   */
  double to_double(const Sdouble &s)
  {
    return s.val();
  }



  /**
   * Return a double given a SSdouble
   */
  double to_double(const SSdouble &s)
  {
    return s.val().val();
  }



  /**
   * Do nothing, required for compatibility.
   */
  template <int index, int dim>
  Tensor <index,dim> to_double(const Tensor<index,dim> &t)
  {
    return t;
  }



  /**
   * Return a Tensor<index,dim,double> given a Tensor<index,dim,Sdouble<number> >
   * where T can be double, Sdouble, etc.
   */
  template <int index, int dim, typename number>
  Tensor <index,dim> to_double(const Tensor<index,dim, Sacado::Fad::DFad<number> > &t)
  {
    Tensor<index,dim> ret;
    for (unsigned int i=0; i<dim; ++i)
      ret[i] = to_double(t[i]);
    return ret;
  }



  /**
   * Return a std::vector<Tensor<index,dim,double> > given a
   * std::vector<Tensor<index,dim,Sdouble<T> > > where T can be
   * double, Sdouble etc.
   */
  template <int index, int dim, typename number>
  std::vector<Tensor<index,dim> > to_double(const std::vector<Tensor<index,dim,Sacado::Fad::DFad<number> > > &v)
  {
    std::vector<Tensor<index,dim> > ret(v.size());
    for (unsigned int q=0; q<v.size(); ++q)
      ret[q] = to_double(v[q]);
    return ret;
  }



  /**
   * Return a std::vector<double> > given a
   * std::vector<number> > > where number can be
   * double, Sdouble etc.
   */
  template <typename number>
  std::vector<double> to_double(const std::vector<number> &v)
  {
    std::vector<double> ret(v.size());
    for (unsigned int q=0; q<v.size(); ++q)
      ret[q] = to_double(v[q]);
    return ret;
  }


}// end namespace


D2K_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS

#endif

