
#include <deal.II/base/function_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iomanip>
#include <vector>

#include "tests.h"
#include "utilities.h"
#include "parameter_acceptor.h"
#include "parsed_dirichlet_bcs.h"


template <int dim>
void test ()
{
  deallog << "dim = " << dim << std::endl;

  Triangulation<dim>   triangulation;
  FE_Q<dim>            fe(1);
  DoFHandler<dim>      dof_handler(triangulation);

  GridGenerator::hyper_cube (triangulation, 0, 1);
  triangulation.begin_active()->face(0)->set_boundary_id(10);
  triangulation.begin_active()->face(1)->set_boundary_id(20);
  triangulation.refine_global (1);

  dof_handler.distribute_dofs (fe);
  QGauss<dim-1> quadrature(2);
  MappingQ1<dim> mapping;

  std::map<types::global_dof_index,double> boundary_values;
//  typename FunctionMap<dim>::type boundary_map;
//  Functions::SquareFunction<dim> f;
//  boundary_map[10] = &f;
//  boundary_map[20] = &f;
//
//  VectorTools::project_boundary_values(mapping, dof_handler, boundary_map, quadrature,
//                                       boundary_values);
  ParsedDirichletBCs<dim,dim,1> parsed_dirichlet("Parsed Dirichlet BCs","","10=0 % 20=0",
                                                 (dim==1 ?"10=x^2 % 20=x^2" :
                                                  (dim==2 ?"10=x^2+y^2 % 20=x^2+y^2" :"10=x^2+y^2+z^2 % 20=x^2+y^2+z^2")));

  ParameterAcceptor::initialize();
  parsed_dirichlet.project_boundary_values(mapping,dof_handler,quadrature,boundary_values);
  deallog << boundary_values.size() << std::endl;
  for (std::map<types::global_dof_index,double>::const_iterator
       p = boundary_values.begin();
       p != boundary_values.end(); ++p)
    deallog << p->first << ' ' << p->second << std::endl;
}


int main ()
{
  initlog();
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  ParameterAcceptor::prm.log_parameters(deallog);

  test<1>();
  test<2>();
  test<3>();
}

