#ifndef _sak_fe_values_cache_h_
#define _sak_fe_values_cache_h_

#include <deal.II/fe/fe_values.h>
#include "sak_data.h"
#include "dof_utilities.h"
#include "utilities.h"


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
    fe_values       (mapping, fe, quadrature, update_flags),
    fe_face_values  (mapping, fe, face_quadrature, face_update_flags),
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
    local_dof_indices(scratch.fe_values.get_fe().dofs_per_cell)
  {};


  /**
   * Initialize the internal FEValues to use the given @p cell.
   */
  void reinit(const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell)
  {
    fe_values.reinit(cell);
    cell->get_dof_indices(local_dof_indices);
  };


  /**
   * Return the quadrature points of the internal FEValues object.
   */
  const std::vector<Point<spacedim> > &
  get_quadrature_points()
  {
    return fe_values.get_quadrature_points();
  }


  /**
   * Return the JxW values of the internal FEValues object.
   */
  const std::vector<double > &
  get_JxW_values()
  {
    return fe_values.get_JxW_values();
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
  void set_solution_vector(const std::string &prefix, const VEC &input_vector, const Number dummy)
  {
    std::string name = prefix+"_independent_local_dofs_"+type(dummy);

    if (!cache.have(name))
      cache.add_copy(std::vector<Number>(fe_values.get_fe().dofs_per_cell), name);

    std::vector<Number> &independent_local_dofs
      = cache.template get<std::vector<Number> >(name);
    DOFUtilities::extract_local_dofs(input_vector, local_dof_indices, independent_local_dofs);
  }




  /**
   * Return the values of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function set_solution_vector with the same name you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <std::vector<Number> > &
  get_values(const std::string &prefix, const Number dummy)
  {

    std::string dofs_name = prefix+"_independent_local_dofs_"+type(dummy);

    std::string name = prefix+"_all_values_"+type(dummy);

    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;
    const unsigned int           n_components  = fe_values.get_fe().n_components();

    Assert(cache.have(dofs_name),
           ExcMessage("You did not call set_solution_vector with the right types!"));

    std::vector<Number> &independent_local_dofs =
      cache.template get<std::vector<Number> >(dofs_name);

    if (!cache.have(name))
      {
        cache.add_copy(std::vector <std::vector<Number> >(n_q_points, std::vector<Number>(n_components)), name);
      }

    std::vector <std::vector<Number> > &ret = cache.template get<std::vector <std::vector<Number> > >(name);
    DOFUtilities::get_values(fe_values, independent_local_dofs, ret);
    return ret;
  }



  /**
   * Return the values of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function set_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Number> &
  get_values(const std::string &prefix,
             const std::string &additional_prefix,
             const FEValuesExtractors::Scalar &scalar_variable,
             const Number dummy)
  {

    std::string dofs_name = prefix+"_independent_local_dofs_"+type(dummy);

    std::string name = prefix+"_"+additional_prefix+"_scalar_values_"+type(dummy);

    const unsigned int           dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int           n_q_points    = fe_values.n_quadrature_points;
    const unsigned int           n_components  = fe_values.get_fe().n_components();

    Assert(cache.have(dofs_name),
           ExcMessage("You did not call set_solution_vector with the right types!"));

    std::vector<Number> &independent_local_dofs =
      cache.template get<std::vector<Number> >(dofs_name);

    if (!cache.have(name))
      {
        cache.add_copy(std::vector<Number>(n_q_points), name);
      }

    std::vector<Number> &ret = cache.template get<std::vector<Number> >(name);
    DOFUtilities::get_values(fe_values, independent_local_dofs, scalar_variable, ret);
    return ret;
  }



  /**
   * Return the values of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function set_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor <1, spacedim, Number> > &
  get_values(const std::string &prefix,
             const std::string &additional_prefix,
             const FEValuesExtractors::Vector &vector_variable,
             const Number dummy)
  {

    // Get the indenpendent dofs
    std::string dofs_name = prefix+"_independent_local_dofs_"+type(dummy);

    Assert(cache.have(dofs_name),
           ExcMessage("You did not call set_solution_vector with the right types!"));

    std::vector<Number> &independent_local_dofs =
      cache.template get<std::vector<Number> >(dofs_name);


    // Now build the return type
    typedef typename std::vector<Tensor<1, spacedim, Number> > RetType;

    std::string name = prefix+"_"+additional_prefix+"_vector_values_"+type(dummy);

    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_values(fe_values, independent_local_dofs, vector_variable, ret);
    return ret;
  }



  /**
   * Return the divergence of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function set_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Number> &
  get_div_values(const std::string &prefix,
                 const std::string &additional_prefix,
                 const FEValuesExtractors::Vector &vector_variable,
                 const Number dummy)
  {

    // Get the indenpendent dofs
    std::string dofs_name = prefix+"_independent_local_dofs_"+type(dummy);

    Assert(cache.have(dofs_name),
           ExcMessage("You did not call set_solution_vector with the right types!"));

    std::vector<Number> &independent_local_dofs =
      cache.template get<std::vector<Number> >(dofs_name);


    // Now build the return type
    typedef typename std::vector<Number> RetType;

    std::string name = prefix+"_"+additional_prefix+"_div_values_"+type(dummy);

    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_div_values(fe_values, independent_local_dofs, vector_variable, ret);
    return ret;
  }



  /**
   * Return the gradient of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function set_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor<1, spacedim, Number> > &
  get_grad_values(const std::string &prefix,
                  const std::string &additional_prefix,
                  const FEValuesExtractors::Scalar &variable,
                  const Number dummy)
  {

    // Get the indenpendent dofs
    std::string dofs_name = prefix+"_independent_local_dofs_"+type(dummy);

    Assert(cache.have(dofs_name),
           ExcMessage("You did not call set_solution_vector with the right types!"));

    std::vector<Number> &independent_local_dofs =
      cache.template get<std::vector<Number> >(dofs_name);


    // Now build the return type
    typedef typename std::vector <Tensor<1, spacedim, Number> > RetType;

    std::string name = prefix+"_"+additional_prefix+"_grad_values_"+type(dummy);

    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_grad_values(fe_values, independent_local_dofs, variable, ret);
    return ret;
  }



  /**
   * Return the gradient of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function set_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor<2, spacedim, Number> > &
  get_grad_values(const std::string &prefix,
                  const std::string &additional_prefix,
                  const FEValuesExtractors::Vector &variable,
                  const Number dummy)
  {

    // Get the indenpendent dofs
    std::string dofs_name = prefix+"_independent_local_dofs_"+type(dummy);

    Assert(cache.have(dofs_name),
           ExcMessage("You did not call set_solution_vector with the right types!"));

    std::vector<Number> &independent_local_dofs =
      cache.template get<std::vector<Number> >(dofs_name);


    // Now build the return type
    typedef typename std::vector <Tensor<2, spacedim, Number> > RetType;

    std::string name = prefix+"_"+additional_prefix+"_grad_values_"+type(dummy);

    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_grad_values(fe_values, independent_local_dofs, variable, ret);
    return ret;
  }


  /**
   * Return the deformation gradient of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function set_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor<2, spacedim, Number> > &
  get_F_values(const std::string &prefix,
               const std::string &additional_prefix,
               const FEValuesExtractors::Vector &variable,
               const Number dummy)
  {

    // Get the indenpendent dofs
    std::string dofs_name = prefix+"_independent_local_dofs_"+type(dummy);

    Assert(cache.have(dofs_name),
           ExcMessage("You did not call set_solution_vector with the right types!"));

    std::vector<Number> &independent_local_dofs =
      cache.template get<std::vector<Number> >(dofs_name);


    // Now build the return type
    typedef typename std::vector <Tensor<2, spacedim, Number> > RetType;

    std::string name = prefix+"_"+additional_prefix+"_F_values_"+type(dummy);

    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_F_values(fe_values, independent_local_dofs, variable, ret);
    return ret;
  }



  /**
   * Return the deformation gradient of the named cached solution vector.
   *
   * Before you call this function, you have to make sure you have previously
   * called the function set_solution_vector with the same @p prefix you give here
   * and with the same dummy type you use here.
   */
  template<typename Number>
  const std::vector <Tensor<2, spacedim, Number> > &
  get_sym_grad_values(const std::string &prefix,
                      const std::string &additional_prefix,
                      const FEValuesExtractors::Vector &variable,
                      const Number dummy)
  {

    // Get the indenpendent dofs
    std::string dofs_name = prefix+"_independent_local_dofs_"+type(dummy);

    Assert(cache.have(dofs_name),
           ExcMessage("You did not call set_solution_vector with the right types!"));

    std::vector<Number> &independent_local_dofs =
      cache.template get<std::vector<Number> >(dofs_name);


    // Now build the return type
    typedef typename std::vector <Tensor<2, spacedim, Number> > RetType;

    std::string name = prefix+"_"+additional_prefix+"_sym_grad_values_"+type(dummy);

    const unsigned int           n_q_points    = fe_values.n_quadrature_points;

    if (!cache.have(name))
      cache.add_copy(RetType(n_q_points), name);

    RetType &ret = cache.template get<RetType>(name);
    DOFUtilities::get_sym_grad_values(fe_values, independent_local_dofs, variable, ret);
    return ret;
  }


private:

  SAKData                                           cache;
  FEValues<dim, spacedim>                           fe_values;
  FEFaceValues<dim, spacedim>                       fe_face_values;
  std::vector<types::global_dof_index>  local_dof_indices;
};

#endif
