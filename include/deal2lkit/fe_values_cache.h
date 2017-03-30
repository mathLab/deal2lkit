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

#ifndef _d2k_fe_values_cache_h_
#define _d2k_fe_values_cache_h_

#include <deal2lkit/config.h>
#include <deal.II/fe/fe_values.h>
#include <deal2lkit/any_data.h>
#include <deal2lkit/dof_utilities.h>
#include <deal2lkit/utilities.h>

D2K_NAMESPACE_OPEN
/**
 * Helper class to simplify the assembly of non-linear problems, and the
 * evaluation of finite element fields. This class will allow you to loop
 * over cells, and access to temporary std::vectors of objects (of the same
 * size of your quadrature formulas, and compatible with the type you specify
 * in the getter functions) without having to declare each of them beforehand,
 * and without having to reinitialize them as Sacado types if the type you
 * require is actually derived by a Sacado type.
 *
 * This is very convenient, for example, if you plan to use Sacado to assemble
 * Jacobian matrices. The following example shows you how to use this class:
 *
 * @code
 * FEValuesCache fev(mapping, fe, quad, flags, face_quad, face_flags);
 *
 * // A couple of extractors, to simplify things.
 * FEValuesExtractors::Vector velocity(0);
 * FEValuesExtractors::Vector pressure(dim);
 *
 * for(auto cell : dh.active_cell_iterators()) {
 *   fev.reinit(cell);
 *
 *   double dummy=0;
 *   fev.cache_local_solution_vector("solution", solution, dummy);
 *   // From now on, you can access to all values/derivatives/divergence of
 *   // the solution finite element field, at the quadrature points
 *
 *
 *   // get the values of the velocity and pressure
 *   // at the quadrature points
 *   std::vector<Tensor<1,dim> > &velocities
 *      = fev.get_values("solution", "v", velocities, dummy);
 *
 *   std::vector<double> &divergences
 *      = fev.get_divergences("solution", "div_v", velocities, dummy);
 *
 *   std::vector<double> &pressures
 *      = fev.get_values("solution", "p", pressure, dummy);
 *
 *   // do something with them
 *   ...
 * }
 * @endcode
 *
 * Internally, FEValuesCache creates unique objects of the right size and type,
 * whenever on of the get_* function is called. These objects are identified by
 * four things: the type of FEValues (either FEValues, or FEFaceValues),
 * the finite element vector identifier ("solution", in the example
 * above), the field type identifier ("v", "div_v", and "p" in the example above),
 * and the type of the dummy variable (used only to determine its type. The value
 * of the variable is ignored).
 *
 * These objects are created the first time these functions are called with a unique
 * identifier, and reused by reference all subsequent times, saving you the necessity
 * to instantiate and initialize new vectors of the right type and size everytime
 * you need to access finite element solutions.
 *
 * Before anything sensible can be extracted by this class, you have to call first
 * the cache_local_solution_vector() function, and then one of the reinit() functions.
 *
 * If you try to access the values before a cache_local_solution_vector is called, an
 * exception will be thrown.
 *
 * This function handles correctly also Sacado types (both first order and second order
 * derivatives), and automatically initializes their derivative. This allows you to
 * access derivatives w.r.t. to
 *
 */
template <int dim, int spacedim=dim>
class FEValuesCache
{
public:
  /**
   * Explicit constructor.
   */
  FEValuesCache     (const Mapping<dim, spacedim>         &mapping,
                     const FiniteElement<dim, spacedim>   &fe,
                     const Quadrature<dim>                &quadrature,
                     const UpdateFlags                    &update_flags,
                     const Quadrature<dim-1>              &face_quadrature,
                     const UpdateFlags                    &face_update_flags):
    fe_values         (mapping, fe, quadrature, update_flags),
    fe_face_values    (mapping, fe, face_quadrature, face_update_flags),
    fe_subface_values (mapping, fe, face_quadrature, face_update_flags),
    local_dof_indices(fe.dofs_per_cell)

  {};

  /**
   * Deep copy constructor.
   */
  FEValuesCache     (const FEValuesCache<dim,spacedim> &scratch):
    cache   (scratch.cache),
    fe_values ( scratch.fe_values.get_mapping(),
                scratch.fe_values.get_fe(),
                scratch.fe_values.get_quadrature(),
                scratch.fe_values.get_update_flags()),
    fe_face_values ( scratch.fe_values.get_mapping(),
                     scratch.fe_values.get_fe(),
                     scratch.fe_face_values.get_quadrature(),
                     scratch.fe_face_values.get_update_flags()),
    fe_subface_values ( scratch.fe_values.get_mapping(),
                        scratch.fe_values.get_fe(),
                        scratch.fe_subface_values.get_quadrature(),
                        scratch.fe_subface_values.get_update_flags()),
    local_dof_indices(scratch.fe_values.get_fe().dofs_per_cell)
  {};


  /**
   * Initialize the internal FEValues to use the given @p cell.
   */
  void reinit(const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell)
  {
    fe_values.reinit(cell);
    cell->get_dof_indices(local_dof_indices);
    cache.template add_ref<FEValuesBase<dim,spacedim> >(fe_values, "FEValuesBase");
  };


  /**
   * Initialize the internal FEFaceValues to use the given @p face_no on the given
   * @p cell.
   */
  void reinit(const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
              const unsigned int face_no)
  {
    fe_face_values.reinit(cell, face_no);
    cell->get_dof_indices(local_dof_indices);
    cache.template add_ref<FEValuesBase<dim,spacedim> >(fe_face_values, "FEValuesBase");
  };


  /**
  * Initialize the internal FESubFaceValues to use the given @p subface_no, on @p face_no,
  * on the given @p cell.
  */
  void reinit(const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
              const unsigned int face_no, const unsigned int subface_no)
  {
    fe_subface_values.reinit(cell, face_no, subface_no);
    cell->get_dof_indices(local_dof_indices);
    cache.template add_ref<FEValuesBase<dim,spacedim> >(fe_subface_values, "FEValuesBase");
  };

  /**
   * @brief Get the reference to @p cache.
   *
   * @return reference to @p cache.
   */
  AnyData &get_cache()
  {
    return cache;
  };

  /**
   * Get the currently initialized FEValues.
   *
   * This function will return the internal FEValues if the reinit(cell) function
   * was called last. If the reinit(cell, face_no) function was called, then this
   * function returns the internal FEFaceValues.
   */
  const FEValuesBase<dim,spacedim> &get_current_fe_values()
  {
    Assert(cache.have("FEValuesBase"),
           ExcMessage("You have to initialize the cache using one of the "
                      "reinit function first!"));
    FEValuesBase<dim,spacedim> &fev =
      cache.template get<FEValuesBase<dim,spacedim> >("FEValuesBase");
    return fev;
  }


  /**
   *
   */
  template<typename Number>
  std::vector<Number> &
  get_current_independent_local_dofs(const std::string &prefix, Number dummy)
  {

    std::string dofs_name = prefix+"_independent_local_dofs_"+type(dummy);

    Assert(cache.have(dofs_name),
           ExcMessage("You did not call cache_local_solution_vector with the right types!"));

    return  cache.template get<std::vector<Number> >(dofs_name);
  }


  /**
   * Return the quadrature points of the internal FEValues object.
   */
  const std::vector<Point<spacedim> > &
  get_quadrature_points()
  {
    return get_current_fe_values().get_quadrature_points();
  }


  /**
   * Return the JxW values of the internal FEValues object.
   */
  const std::vector<double > &
  get_JxW_values()
  {
    return get_current_fe_values().get_JxW_values();
  }


  /**
   * Print the content of the internal cache.
   */
  template <class STREAM>
  void print_info(STREAM &os)
  {
    cache.print_info(os);
  }


  /**
   * Store internally the independent local dof values associated with the
   * internally initialized cell.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function reinit(cell).
   *
   * At every call of this function, a new vector of dof values is generated and
   * stored internally, unless a previous vector with the same name is found.
   * If this is the case, the content of that vector is overwritten.
   *
   * If you give a unique @p prefix, then for each cell you are guaranteed you
   * get a unique vector of independent dofs of the same type as the dummy
   * variable. DOFUtilities is used to initialize the variables, so if you
   * use a Sacado type, this will automatically initialize the derivatives
   * inside the vector.
   */
  template<typename VEC, typename Number>
  void cache_local_solution_vector(const std::string &prefix, const VEC &input_vector, const Number dummy)
  {
    const unsigned int n_dofs = get_current_fe_values().get_fe().dofs_per_cell;

    std::string name = prefix+"_independent_local_dofs_"+type(dummy);

    if (!cache.have(name))
      cache.add_copy(std::vector<Number>(n_dofs), name);

    std::vector<Number> &independent_local_dofs
      = cache.template get<std::vector<Number> >(name);
    DOFUtilities::extract_local_dofs(input_vector, local_dof_indices, independent_local_dofs);
  }




  /**
   * Return the values of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same name you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <std::vector<Number> > &
  get_values(const std::string &prefix, const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;
    const unsigned int           n_components  = fev.get_fe().n_components();

    std::string name = prefix+"_all_values_q"+Utilities::int_to_string(n_q_points)+"_"+type(dummy);

    if (!cache.have(name))
      {
        cache.add_copy(std::vector <std::vector<Number> >(n_q_points, std::vector<Number>(n_components)), name);
      }

    std::vector <std::vector<Number> > &ret = cache.template get<std::vector <std::vector<Number> > >(name);
    DOFUtilities::get_values(fev, independent_local_dofs, ret);
    return ret;
  }



  /**
   * Return the values of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Number> &
  get_values(const std::string &prefix,
             const std::string &additional_prefix,
             const FEValuesExtractors::Scalar &variable,
             const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;

    std::string name = prefix+"_"+additional_prefix+"_scalar_values_q"+
                       Utilities::int_to_string(n_q_points)+"_"+type(dummy);

    // Now build the return type
    typedef typename std::vector<Number> RetType;


    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_values(fev, independent_local_dofs, variable, ret);
    return ret;
  }



  /**
   * Return the values of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor <1, spacedim, Number> > &
  get_values(const std::string &prefix,
             const std::string &additional_prefix,
             const FEValuesExtractors::Vector &vector_variable,
             const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;

    std::string name = prefix+"_"+additional_prefix+"_vector_values_q"+
                       Utilities::int_to_string(n_q_points)+"_"+type(dummy);


    // Now build the return type
    typedef typename std::vector<Tensor<1, spacedim, Number> > RetType;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_values(fev, independent_local_dofs, vector_variable, ret);
    return ret;
  }



  /**
   * Return the divergence of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Number> &
  get_divergences(const std::string &prefix,
                  const std::string &additional_prefix,
                  const FEValuesExtractors::Vector &vector_variable,
                  const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;

    std::string name = prefix+"_"+additional_prefix+"_div_values_q"+
                       Utilities::int_to_string(n_q_points)+"_"+type(dummy);


    // Now build the return type
    typedef typename std::vector<Number> RetType;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_divergences(fev, independent_local_dofs, vector_variable, ret);
    return ret;
  }



  /**
   * Return the gradient of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor<1, spacedim, Number> > &
  get_gradients(const std::string &prefix,
                const std::string &additional_prefix,
                const FEValuesExtractors::Scalar &variable,
                const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;

    std::string name = prefix+"_"+additional_prefix+"_grad_values_q"+
                       Utilities::int_to_string(n_q_points)+"_"+type(dummy);

    // Now build the return type
    typedef typename std::vector <Tensor<1, spacedim, Number> > RetType;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_gradients(fev, independent_local_dofs, variable, ret);
    return ret;
  }



  /**
   * Return the gradient of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor<2, spacedim, Number> > &
  get_gradients(const std::string &prefix,
                const std::string &additional_prefix,
                const FEValuesExtractors::Vector &variable,
                const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;

    std::string name = prefix+"_"+additional_prefix+"_grad2_values_q"+
                       Utilities::int_to_string(n_q_points)+"_"+type(dummy);

    // Now build the return type
    typedef typename std::vector <Tensor<2, spacedim, Number> > RetType;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_gradients(fev, independent_local_dofs, variable, ret);
    return ret;
  }

  /**
   * Return the curl of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor<1, (spacedim > 2 ? spacedim : 1), Number> > &
  get_curls(const std::string &prefix,
            const std::string &additional_prefix,
            const FEValuesExtractors::Vector &variable,
            const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;

    std::string name = prefix+"_"+additional_prefix+"_curl_values_q"+
                       Utilities::int_to_string(n_q_points)+"_"+type(dummy);

    // Now build the return type
    typedef typename std::vector <Tensor<1, (spacedim > 2 ? spacedim : 1), Number> > RetType;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_curls(fev, independent_local_dofs, variable, ret);
    return ret;
  }

  /**
   * Return the laplacian of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Number> &
  get_laplacians(const std::string &prefix,
                 const std::string &additional_prefix,
                 const FEValuesExtractors::Scalar &variable,
                 const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;

    std::string name = prefix+"_"+additional_prefix+"_laplacian_values_q"+
                       Utilities::int_to_string(n_q_points)+"_"+type(dummy);

    // Now build the return type
    typedef typename std::vector <Number> RetType;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_laplacians(fev, independent_local_dofs, variable, ret);
    return ret;
  }


  /**
   * Return the laplacian of the named cached solution vector. Vector value case.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor<1,spacedim,Number> > &
  get_laplacians(const std::string &prefix,
                 const std::string &additional_prefix,
                 const FEValuesExtractors::Vector &variable,
                 const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;

    std::string name = prefix+"_"+additional_prefix+"_laplacian2_values_q"+
                       Utilities::int_to_string(n_q_points)+"_"+type(dummy);

    // Now build the return type
    typedef typename std::vector <Tensor<1,spacedim,Number> > RetType;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_laplacians(fev, independent_local_dofs, variable, ret);
    return ret;
  }

  /**
   * Return the deformation gradient of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor<2, spacedim, Number> > &
  get_deformation_gradients(const std::string &prefix,
                            const std::string &additional_prefix,
                            const FEValuesExtractors::Vector &variable,
                            const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;

    std::string name = prefix+"_"+additional_prefix+"_F_values_q"+
                       Utilities::int_to_string(n_q_points)+"_"+type(dummy);

    // Now build the return type
    typedef typename std::vector <Tensor<2, spacedim, Number> > RetType;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_deformation_gradients(fev, independent_local_dofs, variable, ret);
    return ret;
  }



  /**
   * Return the symmetric gradient of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor<2, spacedim, Number> > &
  get_symmetric_gradients(const std::string &prefix,
                          const std::string &additional_prefix,
                          const FEValuesExtractors::Vector &variable,
                          const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;

    std::string name = prefix+"_"+additional_prefix+"_sym_grad_values_q"+
                       Utilities::int_to_string(n_q_points)+"_"+type(dummy);

    // Now build the return type
    typedef typename std::vector <Tensor<2, spacedim, Number> > RetType;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_symmetric_gradients(fev, independent_local_dofs, variable, ret);
    return ret;
  }

  /**
   * Return the hessian of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor<2,spacedim,Number> > &
  get_hessians(const std::string &prefix,
               const std::string &additional_prefix,
               const FEValuesExtractors::Scalar &variable,
               const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;

    std::string name = prefix+"_"+additional_prefix+"_hessians_values_q"+
                       Utilities::int_to_string(n_q_points)+"_"+type(dummy);

    // Now build the return type
    typedef typename std::vector <Tensor<2,spacedim,Number> > RetType;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_hessians(fev, independent_local_dofs, variable, ret);
    return ret;
  }

  /**
   * Return the hessian of the named cached solution vector. Vector value case.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function cache_local_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor<3,spacedim,Number> > &
  get_hessians(const std::string &prefix,
               const std::string &additional_prefix,
               const FEValuesExtractors::Vector &variable,
               const Number dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_current_independent_local_dofs(prefix, dummy);

    const FEValuesBase<dim,spacedim> &fev = get_current_fe_values();

    const unsigned int           n_q_points    = fev.n_quadrature_points;

    std::string name = prefix+"_"+additional_prefix+"_hessians2_values_q"+
                       Utilities::int_to_string(n_q_points)+"_"+type(dummy);

    // Now build the return type
    typedef typename std::vector <Tensor<3,spacedim,Number> > RetType;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_hessians(fev, independent_local_dofs, variable, ret);
    return ret;
  }

private:

  AnyData                                           cache;
  FEValues<dim, spacedim>                           fe_values;
  FEFaceValues<dim, spacedim>                       fe_face_values;
  FESubfaceValues<dim, spacedim>                    fe_subface_values;
  std::vector<types::global_dof_index>              local_dof_indices;
};

D2K_NAMESPACE_CLOSE

#endif

