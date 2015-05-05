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

#include "error_handler.h"
#include "parsed_grid_generator.h"
#include <deal.II/base/function_lib.h>
#include <deal.II/dofs/dof_handler.h>

int main ()
{
  initlog();

  ErrorHandler<> eh(" ", "u, u, p","L2, Linfty, H1; AddUp; L2"); // Only one table

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

}
