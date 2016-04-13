//-----------------------------------------------------------
//
//    Copyright (C) 2016 by the deal2lkit authors
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

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_math.h>

#include <mpi.h>

#include <nvector/nvector_parallel.h>

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/base/index_set.h>

#include <mpi.h>

#include <deal2lkit/utilities.h>
#include <deal2lkit/kinsol_interface.h>
#include <deal2lkit/parameter_acceptor.h>

#include <fstream>

using namespace deal2lkit;


int main(int argc, char **argv)
{

  Utilities::MPI::MPI_InitFinalize initializater(argc, argv, 1);
  MPI_Comm comm = MPI_COMM_WORLD;

  mpi_initlog();

  int success;
  int numprocs = Utilities::MPI::n_mpi_processes(comm);
  int myid = Utilities::MPI::this_mpi_process(comm);

  typedef TrilinosWrappers::MPI::BlockVector VEC;

  unsigned int ndof = 2;

  // dofs:
  std::function<unsigned int()> n_dofs = [&] ()
  {
    return numprocs*ndof;
  };
  // create new vector:
  std::function<shared_ptr<VEC>()> create_new_vector = [&] ()
  {
    int n_procs = numprocs;
    int my_id = myid;
    static std::vector<IndexSet> all_is;
    if (all_is.size() == 0)
      {
        IndexSet is(n_dofs()/2);
        is.add_range(my_id*(n_dofs()/n_procs/2), (my_id+1)*(n_dofs()/n_procs/2));
        all_is.push_back(is);
        all_is.push_back(is);
      }
    return SP(new TrilinosWrappers::MPI::BlockVector(all_is, comm));
  };
  // residual:
  std::function< int( const VEC &, VEC &) > residual = [&]( const VEC &y, VEC &res )
  {
    res = y;
    return 0;
  };
  // setup_jacobian:
  std::function<int(const VEC &)> setup_jacobian = [](const VEC &)
  {
    return 0;
  };
  // solve_linear_system:
  std::function<int(const VEC &, VEC &)> solve_linear_system = [] (const VEC &, VEC &)
  {
    return 0;
  };
  // Jacobian vector mult:
  std::function<int(const VEC &, VEC &)> jacobian_vmult = [](const VEC &, VEC &)
  {
    return 0;
  };


  // get the Kinsol solver:
  KINSOLInterface<VEC> solver("kinsol", comm );
//  KINSOLInterface<VEC> solver( create_new_vector, residual, setup_jacobian, solve_linear_system, jacobian_vmult, comm );
  // initialize its parameters:
  ParameterAcceptor::initialize(SOURCE_DIR "/parameters/kinsol_interface_trilinos_01.prm", "ode.prm");

  solver.create_new_vector = create_new_vector;
  solver.residual = residual;
  solver.setup_jacobian = setup_jacobian;
  solver.solve_linear_system = solve_linear_system;
  solver.jacobian_vmult = jacobian_vmult;
  // create the solution vector:
  auto solution = solver.create_new_vector();
  // solve the system without changing the strategy and the initial guess:
  success = solver.solve( *solution );
  solution->print(deallog);
  deallog << "On KinSol completion the status flag is: " << success << std::endl;

  // change the residual:
  std::function< int( const VEC &y, VEC &res) > new_residual = [&]( const VEC &y, VEC &res )
  {
    shared_ptr<VEC> temp = create_new_vector();
    *temp = y;
    vector_shift( *temp, 1.0 );
    res = *temp;
    return 0;
  };
  // put the new residual in the solver:
  solver.residual = new_residual;
  // create the solution vector:
  solution = solver.create_new_vector();
  // reinitialize the solver:
  solver.initialize_solver( *solution );
  // solve again:
  success = solver.solve( *solution );

  solution->print(deallog);

  deallog << "On KinSol completion the status flag is: " << success << std::endl;

  return (0);
}
