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
#include <deal.II/base/function.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>

int main ()
{
    initlog();
    
    ParsedGridGenerator<2> gg;
    ErrorHandler<> eh; // Only one table

    ParameterAcceptor::initialize();
    
    auto tria = gg.serial();

    FE_Q<2> fe(1);
    DoFHandler<2> dh(*tria);
    dh.distribute_dofs(fe);
    
}
