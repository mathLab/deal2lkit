#ifndef LH_ERROR_HANDLER_H
#define LH_ERROR_HANDLER_H

#include <fstream>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>

// #include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/config.h>


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe.h>

#include <deal.II/base/parameter_handler.h>

#include "parameter_acceptor.h"
#include <map>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

enum NormFlags
{
  None = 0x00,
  Linfty = 0x01,
  L2 = 0x02,
  W1infty = 0x04,
  H1 = 0x08,
  AddUp = 0x10
};

using namespace dealii;

template <int ntables=1>
class ErrorHandler : public ParameterAcceptor
{
public:
  /** The constructor takes an optional name, specifying the parameter
      entry. */
  ErrorHandler (const std::string name="");

  /** Initialize the given values for the paramter file. */
  virtual void declare_parameters(ParameterHandler &prm);

  /** Parse the given parameter handler. */
  virtual void parse_parameters(ParameterHandler &prm);

  /** Calculate the error of the numeric solution in variuous norms. Store
      the result in the given table. */
  template<int dim, int spacedim, typename VEC>
  void error_from_exact(const DoFHandler<dim,spacedim> &vspace,
                        const VEC &solution,
                        const Function<spacedim> &exact,
                        unsigned int table_no = 0,
                        double dt=0.);

  /** Difference between two solutions in two different vector spaces. */
  template<typename DH, typename VEC>
  void difference(const DH &, const VEC &,
                  const DH &, const VEC &,
                  unsigned int table_no = 0, double dt=0.);

  /** Difference between two solutions in the same vector space. */
  template<typename DH, typename VEC>
  void difference(const DH &, const VEC &,
                  const VEC &, unsigned int table_no = 0, double dt=0.);

  /** By default output first table. */
  void output_table(std::ostream &out=std::cout, const unsigned int table_no=0);

private:
  /** Error results.*/
  std::vector<ConvergenceTable>  tables;

  /** Headers for tables and output. Contains the name of the solution
      components. */
  std::vector<std::string> headers;

  /** Headers for latex tables. Contains the name of the solution
      components. */
  std::vector<std::string> latex_headers;

  /** Captions for latex. */
  std::vector<std::string> latex_captions;

  /** Names of the tables. */
  std::vector<std::string> names;

  /** Type of error to compute per components. */
  std::vector<std::vector<NormFlags> > types;

  /** The parameters have been read. */
  bool initialized;

  /** Compute the error. If this is false, all functions regarding
      errors are disabled and don't do anything.*/
  bool compute_error;

  /** Add convergence rates. */
  std::vector<bool> add_rates;

  /** Write the error files. */
  bool write_error;

  /** Output the error file also on screen. */
  bool output_error;

  /** The error file format. */
  std::string error_file_format;

  /** The extra column to add to the tables. */
  std::vector<std::map<std::string, bool> > extras;

  /** Wether or not to calculate the rates according to the given keys. */
  std::vector<std::string> rate_keys;
};

/**
 * Global operator which returns an object in which all bits are set
 * which are either set in the first or the second argument. This
 * operator exists since if it did not then the result of the bit-or
 * <tt>operator |</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * NormFlags.
 */
inline
NormFlags
operator | (NormFlags f1, NormFlags f2)
{
  return static_cast<NormFlags> (
           static_cast<unsigned int> (f1) |
           static_cast<unsigned int> (f2));
}

/**
 * Global operator which sets the bits from the second argument also
 * in the first one.
 */
inline
NormFlags &
operator |= (NormFlags &f1, NormFlags f2)
{
  f1 = f1 | f2;
  return f1;
}


/**
 * Global operator which returns an object in which all bits are set
 * which are set in the first as well as the second argument. This
 * operator exists since if it did not then the result of the bit-and
 * <tt>operator &</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * NormFlags.
 */
inline
NormFlags
operator & (NormFlags f1, NormFlags f2)
{
  return static_cast<NormFlags> (
           static_cast<unsigned int> (f1) &
           static_cast<unsigned int> (f2));
}


/**
 * Global operator which clears all the bits in the first argument if
 * they are not also set in the second argument.
 */
inline
NormFlags &
operator &= (NormFlags &f1, NormFlags f2)
{
  f1 = f1 & f2;
  return f1;
}

// ============================================================
// Template instantiations
// ============================================================


template <int ntables>
template<typename DH, typename VECTOR>
void ErrorHandler<ntables>::difference(const DH &dh,
                                       const VECTOR &solution1,
                                       const VECTOR &solution2,
                                       unsigned int table_no,
                                       double dt)
{
  AssertThrow(solution1.size() == solution2.size(), ExcDimensionMismatch(
                solution1.size(), solution2.size()));
  VECTOR solution(solution1);
  solution -= solution2;
  error_from_exact(dh, solution,
                   ConstantFunction<DH::dimension>(0, headers.size()), table_no, dt);
}


template <int ntables>
template<int dim, int spacedim, typename VECTOR>
void ErrorHandler<ntables>::error_from_exact(const DoFHandler<dim,spacedim> &dh,
                                             const VECTOR &solution,
                                             const Function<spacedim> &exact,
                                             unsigned int table_no,
                                             double dt)
{
  if (compute_error)
    {
      AssertThrow(initialized, ExcNotInitialized());
      AssertThrow(table_no < types.size(), ExcIndexRange(table_no, 0, names.size()));
      AssertThrow(exact.n_components == types[table_no].size(),
                  ExcDimensionMismatch(exact.n_components, types[table_no].size()));

      std::vector< std::vector<double> > error( exact.n_components, std::vector<double>(4));
      const unsigned int n_active_cells = dh.get_tria().n_active_cells();
      const unsigned int n_dofs=dh.n_dofs();

      if (extras[table_no]["cells"])
        {
          tables[table_no].add_value("cells", n_active_cells);
          tables[table_no].set_tex_caption("cells", "\\# cells");
          tables[table_no].set_tex_format("cells", "r");
        }
      if (extras[table_no]["dofs"])
        {
          tables[table_no].add_value("dofs", n_dofs);
          tables[table_no].set_tex_caption("dofs", "\\# dofs");
          tables[table_no].set_tex_format("dofs", "r");
        }
      if (extras[table_no]["dt"])
        {
          tables[table_no].add_value("dt", dt);
          tables[table_no].set_tex_caption("dt", "\\Delta t");
          tables[table_no].set_tex_format("dt", "r");
        }

      bool compute_Linfty = false;
      bool compute_L2 = false;
      bool compute_W1infty = false;
      bool compute_H1 = false;
      bool add_this = false;

      unsigned int last_non_add = 0;

      for (unsigned int component=0; component < exact.n_components; ++component)
        {
          NormFlags norm = types[table_no][component];

          // Select one Component
          ComponentSelectFunction<dim> select_component ( component, 1. , exact.n_components);

          Vector<float> difference_per_cell (dh.get_tria().n_active_cells());

          QGauss<dim> q_gauss((dh.get_fe().degree+1) * 2);

          // The add bit is set
          add_this = (norm & AddUp);

          if (!add_this)
            {
              last_non_add  = component;
              compute_L2  = ( norm & L2 );
              compute_H1  = ( norm & H1 );
              compute_W1infty = ( norm & W1infty ) ;
              compute_Linfty  = ( norm & Linfty );
            }
          // if add is set, we do not modify the previous selection

          if (compute_L2)
            {
              VectorTools::integrate_difference (//mapping,
                dh, //dof_handler,
                solution,
                exact,
                difference_per_cell,
                q_gauss,
                VectorTools::L2_norm,
                &select_component );
            }

          const double L2_error = difference_per_cell.l2_norm();
          difference_per_cell = 0;

          if (compute_H1)
            {
              VectorTools::integrate_difference (//mapping,
                dh, //dof_handler,
                solution,
                exact,
                difference_per_cell,
                q_gauss,
                VectorTools::H1_norm,
                &select_component );
            }
          const double H1_error = difference_per_cell.l2_norm();
          difference_per_cell = 0;

          if (compute_W1infty)
            {
              VectorTools::integrate_difference (//mapping,
                dh, //dof_handler,
                solution,
                exact,
                difference_per_cell,
                q_gauss,
                VectorTools::W1infty_norm,
                &select_component );
            }

          const double W1inf_error = difference_per_cell.linfty_norm();

          if (compute_Linfty)
            {
              VectorTools::integrate_difference (//mapping,
                dh, //dof_handler,
                solution,
                exact,
                difference_per_cell,
                q_gauss,
                VectorTools::Linfty_norm,
                &select_component );
            }

          const double Linf_error = difference_per_cell.linfty_norm();

          if (add_this)
            {
              AssertThrow(component, ExcMessage("Cannot add on first component!"));

              error[last_non_add][0] = std::max(error[last_non_add][0], Linf_error);
              error[last_non_add][1] += L2_error;
              error[last_non_add][2] = std::max(error[last_non_add][2], W1inf_error);
              error[last_non_add][3] += H1_error;

            }
          else
            {

              error[component][0] = Linf_error;
              error[component][1] = L2_error;
              error[component][2] = W1inf_error;
              error[component][3] = H1_error;

            }
        }

      for (unsigned int j=0; j<exact.n_components; ++j)
        {
          NormFlags norm = types[table_no][j];
          // If this was added, don't do anything
          if ( !(norm & AddUp) )
            {
              if (norm & Linfty)
                {
                  std::string name = headers[j] + "_Linfty";
                  std::string latex_name = "$\\| " +
                                           latex_headers[j] + " - " +
                                           latex_headers[j] +"_h \\|_\\infty $";
                  double this_error =  error[j][0];

                  tables[table_no].add_value(name, this_error);
                  tables[table_no].set_precision(name, 3);
                  tables[table_no].set_scientific(name, true);
                  tables[table_no].set_tex_caption(name, latex_name);
                }

              if (norm & L2)
                {
                  std::string name = headers[j] + "_L2";
                  std::string latex_name = "$\\| " +
                                           latex_headers[j] + " - " +
                                           latex_headers[j] +"_h \\|_0 $";
                  double this_error =  error[j][1];

                  tables[table_no].add_value(name, this_error);
                  tables[table_no].set_precision(name, 3);
                  tables[table_no].set_scientific(name, true);
                  tables[table_no].set_tex_caption(name, latex_name);
                }
              if (norm & W1infty)
                {
                  std::string name = headers[j] + "_W1infty";
                  std::string latex_name = "$\\| " +
                                           latex_headers[j] + " - " +
                                           latex_headers[j] +"_h \\|_{1,\\infty} $";
                  double this_error =  error[j][2];

                  tables[table_no].add_value(name, this_error);
                  tables[table_no].set_precision(name, 3);
                  tables[table_no].set_scientific(name, true);
                  tables[table_no].set_tex_caption(name, latex_name);
                }
              if (norm & H1)
                {
                  std::string name = headers[j] + "_H1";
                  std::string latex_name = "$\\| " +
                                           latex_headers[j] + " - " +
                                           latex_headers[j] +"_h \\|_1 $";
                  double this_error =  error[j][3];

                  tables[table_no].add_value(name, this_error);
                  tables[table_no].set_precision(name, 3);
                  tables[table_no].set_scientific(name, true);
                  tables[table_no].set_tex_caption(name, latex_name);
                }
            }
        }
    }
}



#endif
