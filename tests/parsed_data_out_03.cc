//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------


#include "tests.h"
#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_finite_element.h>
#include <deal2lkit/parsed_grid_generator.h>
#include <deal2lkit/parsed_data_out.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/lac/block_vector.h>

#include <deal.II/base/mpi.h>

template<int dim>
class Test : public ParameterAcceptor
{
public:
  Test();
  virtual void declare_parameters(ParameterHandler &prm);
  void run();

private:

  void make_grid_fe();
  void setup_dofs();
  void output_results();

  ParsedFiniteElement<dim> fe_builder;
  ParsedGridGenerator<dim,dim> tria_builder;
  unsigned int initial_refinement;
  shared_ptr<Triangulation<dim> >      triangulation;
  shared_ptr<FiniteElement<dim,dim> >  fe;
  shared_ptr<DoFHandler<dim> >         dof_handler;
  BlockVector<double> solution;
  ParsedDataOut<dim, dim>                 data_out;
};

template <int dim>
Test<dim>::Test ()
  :
  ParameterAcceptor("Global parameters"),
  tria_builder("Triangulation"),
  fe_builder("FE_Q",
             "FESystem[FE_Q(2)^dim-FE_Q(1)]","u,u,p"),
  data_out("Data out","vtk","","output")
{}

template <int dim>
void Test<dim>::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &initial_refinement, "Initial global refinement", "1",
                Patterns::Integer(0));
}

template <int dim>
void Test<dim>::make_grid_fe ()
{
  ParameterAcceptor::initialize("parameters_ser.prm", "used_parameters_ser.prm");
  triangulation = SP(tria_builder.serial());
  triangulation->refine_global(initial_refinement);
  dof_handler = SP(new DoFHandler<dim>(*triangulation));
  fe=SP(fe_builder());
}

template <int dim>
void Test<dim>::setup_dofs ()
{
  std::vector<unsigned int> sub_blocks (dim+1,0);
  sub_blocks[dim] = 1;
  dof_handler->distribute_dofs (*fe);

  std::vector<types::global_dof_index> dofs_per_block (2);
  DoFTools::count_dofs_per_block (*dof_handler, dofs_per_block,
                                  sub_blocks);

  const unsigned int n_u = dofs_per_block[0],
                     n_p = dofs_per_block[1];

  solution.reinit (2);
  solution.block(0).reinit (n_u);
  solution.block(1).reinit (n_p);
  solution.collect_sizes ();

  std::vector<unsigned int> block_component (dim+1,0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise (*dof_handler, block_component);

  solution.block(0) = 1.;
}

template <int dim>
void Test<dim>::output_results()
{
  data_out.prepare_data_output( *dof_handler);
  data_out.add_data_vector(solution, fe_builder.get_component_names());
  data_out.write_data_and_clear();
  append_to_file("output.vtk","output");
}

template <int dim>
void Test<dim>::run()
{
  make_grid_fe();
  setup_dofs();
  output_results();
}

int main (int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);
  mpi_initlog();
#else
  initlog();
#endif

  const int dim = 2;

  Test<dim> test;
  test.run();

  // return 0;
}
