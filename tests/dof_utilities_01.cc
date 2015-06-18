// test the DOFUtilities functions
// for Number=double

#include "tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <fstream>

#include "dof_utilities.h"


template<int dim>
void test (const Triangulation<dim> &tr,
           const FiniteElement<dim> &fe)
{
  deallog << "FE=" << fe.get_name()
          << std::endl;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);


  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (fe, quadrature,
                           update_values | update_gradients);

  std::vector<types::global_dof_index>    local_dof_indices (fe_values.dofs_per_cell);
  std::vector<double> independent_local_dof_values (fe_values.dofs_per_cell);

  fe_values.reinit (dof.begin_active());
  dof.begin_active()->get_dof_indices (local_dof_indices);

  Vector<double> global_vector(dof.n_dofs());
  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    global_vector[i] += i*i;

  DOFUtilities::extract_local_dofs(global_vector, local_dof_indices, independent_local_dof_values);

  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    deallog << independent_local_dof_values[i] << std::endl;



  std::vector <double> scalar_values(quadrature.size());
  std::vector <Tensor <1, dim, double> > grad_s(quadrature.size());

  std::vector <double> div_values(quadrature.size());
  std::vector <Tensor <1, dim, double> > vector_values(quadrature.size());
  std::vector <Tensor <2, dim, double> > grad_v(quadrature.size());
  std::vector <Tensor <2, dim, double> > sym_grad_v(quadrature.size());
  FEValuesExtractors::Scalar scalar (dim);
  FEValuesExtractors::Vector vector (0);

  DOFUtilities::get_values(fe_values, independent_local_dof_values, scalar, scalar_values);
  DOFUtilities::get_grad_values(fe_values, independent_local_dof_values, scalar, grad_s);

  DOFUtilities::get_values(fe_values, independent_local_dof_values, vector, vector_values);
  DOFUtilities::get_div_values(fe_values, independent_local_dof_values, vector, div_values);
  DOFUtilities::get_grad_values(fe_values, independent_local_dof_values, vector, grad_v);
  DOFUtilities::get_sym_grad_values(fe_values, independent_local_dof_values, vector, sym_grad_v);

  deallog << "scalar_values" << std::endl;
  for (unsigned int q=0; q<quadrature.size(); ++q)
    deallog << scalar_values[q] << std::endl;

  deallog <<std::endl;
  deallog << "div_values" << std::endl;
  for (unsigned int q=0; q<quadrature.size(); ++q)
    deallog << div_values[q] << std::endl;

  deallog <<std::endl;
  deallog << "grad of a scalar field" << std::endl;
  for (unsigned int q=0; q<quadrature.size(); ++q)
    {
      for (unsigned int d=0; d<dim; ++d)
        deallog << grad_s[q][d] << " ";
      deallog <<std::endl;
    }

  deallog <<std::endl;
  deallog << "vector_values" << std::endl;
  for (unsigned int q=0; q<quadrature.size(); ++q)
    {
      for (unsigned int d=0; d<dim; ++d)
        deallog << vector_values[q][d] << " ";
      deallog <<std::endl;
    }

  deallog <<std::endl;
  deallog << "grad of a vector field" << std::endl;
  for (unsigned int q=0; q<quadrature.size(); ++q)
    {
      for (unsigned int d=0; d<dim; ++d)
        {
          for (unsigned int dd=0; dd<dim; ++dd)
            deallog << grad_v[q][d][d] << " ";
          deallog <<std::endl;
        }
      deallog <<std::endl;
    }

  deallog <<std::endl;
  deallog << "sym grad of a vector field" << std::endl;
  for (unsigned int q=0; q<quadrature.size(); ++q)
    {
      for (unsigned int d=0; d<dim; ++d)
        {
          for (unsigned int dd=0; dd<dim; ++dd)
            deallog << sym_grad_v[q][d][d] << " ";
          deallog <<std::endl;
        }
      deallog <<std::endl;
    }

}



template<int dim>
void test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);


  FESystem<dim> fe (FE_Q<dim>(2), dim,
                    FE_Q<dim>(1), 1);
  test(tr, fe);
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test_hyper_cube<2>();
  test_hyper_cube<3>();
}
