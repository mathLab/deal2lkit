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

#include <deal.II/base/mpi.h>
#include "utilities.h"
#include "parameter_acceptor.h"
#include <iostream>

using namespace dealii;

class ProblemParameters : public ParameterAcceptor
{

  virtual void declare_parameters(ParameterHandler &prm)
  {
    add_parameter(prm, &n_refinements, "Refinement levels", "5", Patterns::Integer());
  };

  int n_refinements;
};


int main (int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize (argc, argv, numbers::invalid_unsigned_int);

      ParameterHandler prm;
      ProblemParameters parameter_class;
      ParameterAcceptor::declare_all_parameters(prm);
      prm.read_input("template_parameters.prm");
      ParameterAcceptor::parse_all_parameters(prm);

      // Do something with your class
      // ProblemClass pb;
      // pb.run();

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
