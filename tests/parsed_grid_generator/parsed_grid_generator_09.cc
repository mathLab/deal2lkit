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


#include "../tests.h"
#include <deal2lkit/utilities.h>
#include <deal2lkit/parsed_grid_generator.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/base/utilities.h>


using namespace deal2lkit;

template<int dim, int spacedim>
void test(ParsedGridGenerator<dim, spacedim> &pgg)
{
  Triangulation<dim, spacedim> *tria = pgg.serial();
  GridOut go;
  go.write_msh(*tria, deallog.get_file_stream());
  delete tria;
}


int main ()
{
  initlog();
  ParsedGridGenerator<2,2> a("Unit Hyperball", "hyper_ball");
  ParsedGridGenerator<2,2> b("Sub Hyper Rectangle", "rectangle");
  ParsedGridGenerator<2,2> c("Hyper Shell", "hyper_shell", "", "", "", "false", "1.0", "0.5", "1.5",
                             "0", "0","","none","");
  ParsedGridGenerator<1,2> d("Hyper Sphere", "hyper_sphere");
  ParsedGridGenerator<3,3> e("Hyper Cube with Cylindrical Hole", "hyper_cube_with_cylindrical_hole");
  ParsedGridGenerator<3,3> f("Hyper L", "hyper_L");
  ParsedGridGenerator<3,3> g("Half Hyper Ball", "half_hyper_ball");
  ParsedGridGenerator<2,2> h("Cylinder", "cylinder");
  ParsedGridGenerator<2,2> i("Truncated Cone", "truncated_cone");
  // ParsedGridGenerator<2,2> l("Hyper Cross", "hyper_cross");
  ParsedGridGenerator<2,2> m("Hyper Cube Slit", "hyper_cube_slit");
  ParsedGridGenerator<2,2> n("Half Hyper Shell", "half_hyper_shell");
  ParsedGridGenerator<2,2> o("Quarter Hyper Shell", "quarter_hyper_shell");
  ParsedGridGenerator<3,3> p("Cylinder Shell", "cylinder_shell", "", "", "", "false", "1.0", "0.5", "1.5",
                             "0", "0","","none","");
  ParsedGridGenerator<2,3> q("Torus", "torus");
  // ParsedGridGenerator<3,3> q("Moebius", "moebius");
  ParsedGridGenerator<2,2> s("Cheese", "cheese", "", "", "", "false", "1.0", "0.5", "1.5",
                             "0", "0","1,2","none","");

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  deallog <<"Unit Hyperball"<<std::endl;
  test(a);
  deallog <<"Sub Hyper Rectangle"<<std::endl;
  test(b);
  deallog <<"Hyper Shell"<<std::endl;
  test(c);
  deallog <<"Hyper Sphere"<<std::endl;
  test(d);
  deallog <<"Hyper L"<<std::endl;
  test(f);
  deallog <<"Half Hyper Ball"<<std::endl;
  test(g);
  deallog <<"Cylinder"<<std::endl;
  test(h);
  deallog <<"Truncated Cone"<<std::endl;
  test(i);
  // deallog <<"Hyper Cross"<<std::endl;
  // test(l);
  deallog <<"Hyper Cube Slit"<<std::endl;
  test(m);
  deallog <<"Half Hyper Shell"<<std::endl;
  test(n);
  deallog <<"Quarter Hyper Shell"<<std::endl;
  test(o);
  // deallog <<"Cylinder Shell"<<std::endl;
  // test(p); // Debug and Release mode have different outputs.
  deallog <<"Torus"<<std::endl;
  test(q);
  // deallog <<"Moebius"<<std::endl;
  // test(r);
  deallog <<"Cheese"<<std::endl;
  test(s);
}
