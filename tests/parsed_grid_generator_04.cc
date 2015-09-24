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
  ParsedGridGenerator<2,2> a("Rectangle");
  ParsedGridGenerator<2,2> b("Rectangle with points");

  ParameterHandler prm;

  Point<2> p1, p2;
  std::string string1, string2, string3;

  p1[0] = -1.;
  p1[1] = -2.;
  p2[0] = 1.;
  p2[1] = 2.;


  string1 = b.create_default_value(p1);
  string2 = b.create_default_value(p2);
  ParameterAcceptor::declare_all_parameters(prm);

  prm.read_input_from_string(""
                             "subsection Rectangle\n"
                             "  set Grid to generate = rectangle \n"
                             "  set Optional Point<spacedim> 1 = 0, 0 \n"
                             "  set Optional Point<spacedim> 2 = -1., -2. \n"
                             "end\n");

  string3 = ""
            "subsection Rectangle with points\n"
            "  set Grid to generate = hyperrectangle \n"
            "  set Optional Point<spacedim> 1 = "+ string1 +"\n"
            "  set Optional Point<spacedim> 2 = "+ string2 +"\n"
            "end\n";
  prm.read_input_from_string(string3.c_str());

  ParameterAcceptor::parse_all_parameters(prm);

  test(a);
  test(b);
}
