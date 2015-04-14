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
#include "utilities.h"
#include "parsed_grid_generator.h"

#include <deal.II/grid/grid_out.h>
#include <deal.II/base/utilities.h>

template<int dim, int spacedim>
void test(ParsedGridGenerator<dim, spacedim> &pgg)
{
  auto tria = pgg.serial();
  GridOut go;
  go.write_msh(*tria, deallog.get_file_stream());
  std::ofstream ofile(("/tmp/mesh_"+Utilities::int_to_string(dim)
                       +Utilities::int_to_string(spacedim)+".msh").c_str());
}

int main ()
{
  initlog();
  ParsedGridGenerator<3,3> a("Subdivided HyperCube");
  ParsedGridGenerator<2,3> b("Subdivided HyperCube codim 1");

  ParameterHandler prm;

  std::vector<unsigned int> int_vec_2(2);
  std::vector<unsigned int> int_vec_3(3);
  std::string int_string_2, int_string_3, string3, string2;

  int_vec_2[0] = 2;
  int_vec_2[1] = 2;
  int_vec_3[0] = 1;
  int_vec_3[1] = 2;
  int_vec_3[2] = 2;


  int_string_2 = b.create_default_value(int_vec_2);
  int_string_3 = b.create_default_value(int_vec_3);
  ParameterAcceptor::declare_all_parameters(prm);

  string3 = ""
    "subsection Subdivided HyperCube\n"
    "  set Grid to generate = rectangle\n"
    "  set Optional Point<spacedim> 1 = -1.,-1., -1.\n"
    "  set Optional Point<spacedim> 2 =  1., 1.,  1.\n"
    "  set Optional vector of dim int = "+ int_string_3 + "\n" 
    "end\n";

  string2 = ""
    "subsection Subdivided HyperCube codim 1\n"
    "  set Grid to generate = rectangle\n"
    "  set Optional Point<spacedim> 1 = -1.,-1., 0.\n"
    "  set Optional Point<spacedim> 2 =  1., 1., 0.\n"
    "  set Optional vector of dim int = "+ int_string_2 + "\n" 
    "end\n";
  
  prm.read_input_from_string(string3.c_str());
  prm.read_input_from_string(string2.c_str());

  ParameterAcceptor::parse_all_parameters(prm);

  test(a);
  test(b);
}
