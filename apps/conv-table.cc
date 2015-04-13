/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth and Ralf Hartmann, University of Heidelberg, 2000
 */



#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/fe/fe_values.h>

#include <typeinfo>
#include <fstream>
#include <iostream>

#include "parsed_grid_generator.h"
#include "parsed_finite_element.h"
#include "error_handler.h"
#include "utilities.h"

namespace ConvTables
{
  using namespace dealii;


  template <int dim>
  class SolutionBase
  {
  protected:
    static const unsigned int n_source_centers = 3;
    static const Point<dim>   source_centers[n_source_centers];
    static const double       width;
  };


  template <>
  const Point<1>
  SolutionBase<1>::source_centers[SolutionBase<1>::n_source_centers]
    = { Point<1>(-1.0 / 3.0),
        Point<1>(0.0),
        Point<1>(+1.0 / 3.0)
      };

  template <>
  const Point<2>
  SolutionBase<2>::source_centers[SolutionBase<2>::n_source_centers]
    = { Point<2>(-0.5, +0.5),
        Point<2>(-0.5, -0.5),
        Point<2>(+0.5, -0.5)
      };

  template <int dim>
  const double SolutionBase<dim>::width = 1./8.;



  template <int dim>
  class Solution : public Function<dim>,
    protected SolutionBase<dim>
  {
  public:
    Solution () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;
  };


  template <int dim>
  double Solution<dim>::value (const Point<dim>   &p,
                               const unsigned int) const
  {
    double return_value = 0;
    for (unsigned int i=0; i<this->n_source_centers; ++i)
      {
        const Tensor<1,dim> x_minus_xi = p - this->source_centers[i];
        return_value += std::exp(-x_minus_xi.norm_square() /
                                 (this->width * this->width));
      }

    return return_value;
  }


  template <int dim>
  Tensor<1,dim> Solution<dim>::gradient (const Point<dim>   &p,
                                         const unsigned int) const
  {
    Tensor<1,dim> return_value;

    for (unsigned int i=0; i<this->n_source_centers; ++i)
      {
        const Tensor<1,dim> x_minus_xi = p - this->source_centers[i];

        return_value += (-2 / (this->width * this->width) *
                         std::exp(-x_minus_xi.norm_square() /
                                  (this->width * this->width)) *
                         x_minus_xi);
      }

    return return_value;
  }



  template <int dim>
  class RightHandSide : public Function<dim>,
    protected SolutionBase<dim>
  {
  public:
    RightHandSide () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
  };


  template <int dim>
  double RightHandSide<dim>::value (const Point<dim>   &p,
                                    const unsigned int) const
  {
    double return_value = 0;
    for (unsigned int i=0; i<this->n_source_centers; ++i)
      {
        const Tensor<1,dim> x_minus_xi = p - this->source_centers[i];

        return_value += ((2*dim - 4*x_minus_xi.norm_square()/
                          (this->width * this->width)) /
                         (this->width * this->width) *
                         std::exp(-x_minus_xi.norm_square() /
                                  (this->width * this->width)));
        return_value += std::exp(-x_minus_xi.norm_square() /
                                 (this->width * this->width));
      }

    return return_value;
  }



  template <int dim>
  class HelmholtzProblem
  {
  public:
    enum RefinementMode
    {
      global_refinement=0, adaptive_refinement=1
    };

    HelmholtzProblem (const RefinementMode      refinement_mode);

    ~HelmholtzProblem ();

    void run ();

  private:
    void make_grid_fe();
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void process_solution (const unsigned int cycle);

    SmartPointer<Triangulation<dim> >    triangulation;
    SmartPointer<FiniteElement<dim,dim> >             fe;
    SmartPointer<DoFHandler<dim> >       dof_handler;


    ConstraintMatrix                        hanging_node_constraints;

    SparsityPattern                         sparsity_pattern;
    SparseMatrix<double>                    system_matrix;

    Vector<double>                          solution;
    Vector<double>                          system_rhs;

    const RefinementMode                    refinement_mode;

    // ConvergenceTable                        convergence_table;

    ErrorHandler<2> eh;
    ParsedGridGenerator<dim,dim> pgg;
    ParsedFiniteElement<dim,dim> fe_builder;
  };




  template <int dim>
  HelmholtzProblem<dim>::HelmholtzProblem (const RefinementMode refinement_mode) :
    refinement_mode (refinement_mode),
//    eh("Error Tables"),
    pgg("Grid"),
    fe_builder("Finite Element")
  {}



  template <int dim>
  HelmholtzProblem<dim>::~HelmholtzProblem ()
  {
    smart_delete(dof_handler);
    smart_delete(fe);
    smart_delete(triangulation);
  }

  template <int dim>
  void HelmholtzProblem<dim>::make_grid_fe ()
  {


    ParameterAcceptor::initialize("params_conv.prm");

    triangulation = pgg.serial();
    dof_handler = new DoFHandler<dim>(*triangulation);
    //GridGenerator::hyper_cube (triangulation, -1, 1);

    std::cout << "Number of active cells: "
              << triangulation->n_active_cells()
              << std::endl;

    fe=fe_builder();

  }


  template <int dim>
  void HelmholtzProblem<dim>::setup_system ()
  {
    dof_handler->distribute_dofs (*fe);
    DoFRenumbering::Cuthill_McKee (*dof_handler);

    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (*dof_handler,
                                             hanging_node_constraints);
    hanging_node_constraints.close ();

    DynamicSparsityPattern dsp (dof_handler->n_dofs(), dof_handler->n_dofs());
    DoFTools::make_sparsity_pattern (*dof_handler, dsp);
    hanging_node_constraints.condense (dsp);
    sparsity_pattern.copy_from (dsp);

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dof_handler->n_dofs());
    system_rhs.reinit (dof_handler->n_dofs());
  }



  template <int dim>
  void HelmholtzProblem<dim>::assemble_system ()
  {
    QGauss<dim>   quadrature_formula(3);
    QGauss<dim-1> face_quadrature_formula(3);

    const unsigned int n_q_points    = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    const unsigned int dofs_per_cell = fe->dofs_per_cell;

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    FEValues<dim>  fe_values (*fe, quadrature_formula,
                              update_values   | update_gradients |
                              update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face_values (*fe, face_quadrature_formula,
                                      update_values         | update_quadrature_points  |
                                      update_normal_vectors | update_JxW_values);

    const RightHandSide<dim> right_hand_side;
    std::vector<double>  rhs_values (n_q_points);

    const Solution<dim> exact_solution;

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler->begin_active(),
    endc = dof_handler->end();
    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit (cell);

        right_hand_side.value_list (fe_values.get_quadrature_points(),
                                    rhs_values);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += ((fe_values.shape_grad(i,q_point) *
                                      fe_values.shape_grad(j,q_point)
                                      +
                                      fe_values.shape_value(i,q_point) *
                                      fe_values.shape_value(j,q_point)) *
                                     fe_values.JxW(q_point));

              cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                              rhs_values [q_point] *
                              fe_values.JxW(q_point));
            }

        for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
          if (cell->face(face_number)->at_boundary()
              &&
              (cell->face(face_number)->boundary_id() == 1))
            {
              fe_face_values.reinit (cell, face_number);

              for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
                {
                  const double neumann_value
                    = (exact_solution.gradient (fe_face_values.quadrature_point(q_point)) *
                       fe_face_values.normal_vector(q_point));

                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    cell_rhs(i) += (neumann_value *
                                    fe_face_values.shape_value(i,q_point) *
                                    fe_face_values.JxW(q_point));
                }
            }

        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              system_matrix.add (local_dof_indices[i],
                                 local_dof_indices[j],
                                 cell_matrix(i,j));

            system_rhs(local_dof_indices[i]) += cell_rhs(i);
          }
      }

    hanging_node_constraints.condense (system_matrix);
    hanging_node_constraints.condense (system_rhs);

    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (*dof_handler,
                                              0,
                                              Solution<dim>(),
                                              boundary_values);
    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);
  }



  template <int dim>
  void HelmholtzProblem<dim>::solve ()
  {
    SolverControl           solver_control (1000, 1e-12);
    SolverCG<>              cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve (system_matrix, solution, system_rhs,
              preconditioner);

    hanging_node_constraints.distribute (solution);
  }



  template <int dim>
  void HelmholtzProblem<dim>::refine_grid ()
  {
    switch (refinement_mode)
      {
      case global_refinement:
      {
        triangulation->refine_global (1);
        break;
      }

      case adaptive_refinement:
      {
        Vector<float> estimated_error_per_cell (triangulation->n_active_cells());

        KellyErrorEstimator<dim>::estimate (*dof_handler,
                                            QGauss<dim-1>(3),
                                            typename FunctionMap<dim>::type(),
                                            solution,
                                            estimated_error_per_cell);

        GridRefinement::refine_and_coarsen_fixed_number (*triangulation,
                                                         estimated_error_per_cell,
                                                         0.3, 0.03);

        triangulation->execute_coarsening_and_refinement ();

        break;
      }

      default:
      {
        Assert (false, ExcNotImplemented());
      }
      }
  }



  template <int dim>
  void HelmholtzProblem<dim>::process_solution (const unsigned int cycle)
  {
    eh.error_from_exact(*dof_handler, solution, Solution<dim>(), refinement_mode);
  }



  template <int dim>
  void HelmholtzProblem<dim>::run ()
  {
    const unsigned int n_cycles = (refinement_mode==global_refinement)?5:9;
    for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
      {
        if (cycle == 0)
          {
            make_grid_fe();
            triangulation->refine_global (3);

            typename Triangulation<dim>::cell_iterator
            cell = triangulation->begin (),
            endc = triangulation->end();
            for (; cell!=endc; ++cell)
              for (unsigned int face_number=0;
                   face_number<GeometryInfo<dim>::faces_per_cell;
                   ++face_number)
                if ((std::fabs(cell->face(face_number)->center()(0) - (-1)) < 1e-12)
                    ||
                    (std::fabs(cell->face(face_number)->center()(1) - (-1)) < 1e-12))
                  cell->face(face_number)->set_boundary_id (1);
          }
        else
          refine_grid ();


        setup_system ();

        assemble_system ();
        solve ();

        process_solution (cycle);
      }

    eh.output_table(std::cout, refinement_mode);

    std::string vtk_filename;
    switch (refinement_mode)
      {
      case global_refinement:
        vtk_filename = "solution-global";
        break;
      case adaptive_refinement:
        vtk_filename = "solution-adaptive";
        break;
      default:
        Assert (false, ExcNotImplemented());
      }

    switch (fe->degree)
      {
      case 1:
        vtk_filename += "-q1";
        break;
      case 2:
        vtk_filename += "-q2";
        break;

      default:
        Assert (false, ExcNotImplemented());
      }

    vtk_filename += ".vtk";
    std::ofstream output (vtk_filename.c_str());

    DataOut<dim> data_out;
    data_out.attach_dof_handler (*dof_handler);
    data_out.add_data_vector (solution, "solution");

    data_out.build_patches (fe->degree);
    data_out.write_vtk (output);


  }

}


int main ()
{
  const unsigned int dim = 2;

  try
    {
      using namespace dealii;
      using namespace ConvTables;

      deallog.depth_console (0);

      {
        HelmholtzProblem<dim>
        helmholtz_problem_2d (HelmholtzProblem<dim>::adaptive_refinement);

        helmholtz_problem_2d.run ();

        std::cout << std::endl;
      }
      ParameterAcceptor::clear();
      {
        HelmholtzProblem<dim>
        helmholtz_problem_2d (HelmholtzProblem<dim>::global_refinement);

        helmholtz_problem_2d.run ();

        std::cout << std::endl;
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}


namespace ConvTables
{
  template const double SolutionBase<2>::width;
}
