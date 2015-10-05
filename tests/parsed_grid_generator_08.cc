//-----------------------------------------------------------
//
//    Copyright (C) 2015 by the deal2lkit authors
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

#include <deal.II/grid/grid_out.h>
#include <deal.II/base/utilities.h>

template<int dim, int spacedim>
void test(ParsedGridGenerator<dim, spacedim> &pgg)
{
  Triangulation<dim, spacedim> *tria = pgg.serial();
  GridOut go;
  go.write_msh(*tria, deallog.get_file_stream());
}


int main ()
{
  initlog();
  ParsedGridGenerator<2,2> a("Rectangle", "rectangle","","10,10", "20,20", "true");
  ParsedGridGenerator<3,3> b("Cube", "rectangle", "", "7,8,9", "15,16,16");


  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  deallog <<"2D"<<std::endl;
  test(a);
  deallog <<"3D"<<std::endl;
  test(b);
}

