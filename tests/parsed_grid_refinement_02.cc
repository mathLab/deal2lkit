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
#include <deal2lkit/parsed_grid_generator.h>
#include <deal2lkit/parsed_grid_refinement.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>

template<int dim, int spacedim>
void test()
{
  ParsedGridGenerator<dim, spacedim> pgg;
  ParsedGridRefinement pgr("", "number");

  ParameterAcceptor::initialize();

  Triangulation<dim, spacedim> *tria = pgg.serial();

  tria->refine_global(3);

  Vector<float> criteria(tria->n_active_cells());

  for (auto cell : tria->active_cell_iterators())
    criteria[cell->index()] = cell->center().norm();

  pgr.mark_cells(criteria, *tria);

  tria->execute_coarsening_and_refinement();

  GridOut go;
  go.write_msh(*tria, deallog.get_file_stream());
  std::ofstream ofile(("/tmp/mesh_"+Utilities::int_to_string(dim)
                       +Utilities::int_to_string(spacedim)+".msh").c_str());
  go.write_msh(*tria, ofile);

  delete tria;
}

int main ()
{
  initlog();

  test<1,1>();
  test<1,2>();
  test<2,2>();
  test<2,3>();
  test<3,3>();
}
