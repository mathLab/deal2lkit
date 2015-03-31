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
void test(ParsedGridGenerator<dim, spacedim> &pgg) {
  auto tria = pgg.serial();
  GridOut go;
  go.write_msh(*tria, deallog.get_file_stream());
  std::ofstream ofile(("/tmp/mesh_"+Utilities::int_to_string(dim)
		       +Utilities::int_to_string(spacedim)+".msh").c_str());
}

int main ()
{
    initlog();
    ParsedGridGenerator<2,2> a("Unit cube");
    ParsedGridGenerator<2,2> b();
    
    ParameterHandler prm;
    ParameterAcceptor::declare_all_parameters(prm);
    prm.read_input_from_string(""
                               "subsection ParsedGridGenerator<2,2>\n"
			       "  set Grid name = file \n"
                               "  set Grid file name = 1.0\n"
                               "end\n");

    prm.log_parameters(deallog);
    ParameterAcceptor::parse_all_parameters(prm);

    test(a);
}
