#include "tests.h"

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

#include "sundials_interface.h"
#include "dae_time_integrator.h"
#include "parameter_acceptor.h"

#include <fstream>

template <typename VEC>
class Solver : public SundialsInterface<VEC>, public ParameterAcceptor
{
public:

  Solver(const MPI_Comm &comm) :
    SundialsInterface<VEC>(comm),
    n_dofs_(1),
    ofile("output.gpl")
  {};


  virtual void declare_parameters(ParameterHandler &prm)
  {
    add_parameter(prm, &n_dofs_, "Number of comps", "2", Patterns::Integer(0));
  }

  virtual unsigned int n_dofs() const
  {
    static int n_procs = Utilities::MPI::n_mpi_processes(this->get_comm());
    static int my_id = Utilities::MPI::this_mpi_process(this->get_comm());
    return n_procs*n_dofs_;
  };

  virtual shared_ptr<VEC> create_new_vector() const
  {
    static int n_procs = Utilities::MPI::n_mpi_processes(this->get_comm());
    static int my_id = Utilities::MPI::this_mpi_process(this->get_comm());
    static std::vector<IndexSet> all_is;
    if (all_is.size() == 0)
      {
        IndexSet is(n_dofs()/2);
        is.add_range(my_id*(n_dofs()/n_procs/2), (my_id+1)*(n_dofs()/n_procs/2));
        all_is.push_back(is);
        all_is.push_back(is);
      }
    return SP(new TrilinosWrappers::MPI::BlockVector(all_is, this->get_comm()));
  }

  virtual void output_step(const double t,
                           const VEC &solution,
                           const VEC &solution_dot,
                           const unsigned int step_number,
                           const double h)
  {
    deallog << "Step " << step_number
            << ", t: " << t << std::endl;
    if (solution.block(0).in_local_range(0))
      {
        ofile << t << " " << solution.block(0)[0]
              << " " << solution.block(1)[0] << std::endl;
      }
  }


  virtual int residual(const double t,
                       const VEC &src_yy,
                       const VEC &src_yp,
                       VEC &dst)
  {

    auto &y = src_yy.block(0);
    auto &p = src_yy.block(1);
    auto &ydot = src_yp.block(0);

    auto &dy = dst.block(0);
    auto &dp = dst.block(1);

    // y_dot - y * p
    dy = y;
    dy.scale(p);
    dy += ydot;

    dp = p;
    dp.add(-1);
    dp.scale(p);
    return 0;
  }

  /** Jacobian vector product. */
  virtual int setup_jacobian(const double t,
                             const VEC &src_yy,
                             const VEC &src_yp,
                             const VEC &residual,
                             const double alpha_in)
  {
    y = SP(new VEC(src_yy));
    y_dot = SP(new VEC(src_yp));
    alpha = alpha_in;
    return 0;
  }


  virtual VEC &differential_components() const
  {
    static shared_ptr<VEC> diff = create_new_vector();
    static bool initialized = false;
    if (initialized == false)
      {
        diff->block(0) = 1.0;
        initialized = true;
      }
    return *diff;
  }

  virtual VEC &get_local_tolerances() const
  {
    static shared_ptr<VEC> diff = create_new_vector();
    static bool initialized = false;
    if (initialized == false)
      {
        diff->add(1.0);
        initialized = true;
      }
    return *diff;
  }

  virtual int solve_jacobian_system(const double t,
                                    const VEC &y,
                                    const VEC &y_dot,
                                    const VEC &residual,
                                    const double alpha,
                                    const VEC &src,
                                    VEC &dst) const
  {
    deallog<<"solved" << std::endl;
    return 0;
  }


private:
  unsigned int n_dofs_;
  std::ofstream ofile;

  shared_ptr<VEC> y;
  shared_ptr<VEC> y_dot;
  double alpha;

};



int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize initializater(argc, argv, 1);

  MPI_Comm comm = MPI_COMM_WORLD;

  mpi_initlog();

  int numprocs = Utilities::MPI::n_mpi_processes(comm);
  int myid = Utilities::MPI::this_mpi_process(comm);
//
//  std::cout << "Process " << getpid() << " is " << myid
//            << " of " << numprocs << " processes" << std::endl;
//  if (myid == 0) system("read -p \"Press [Enter] key to start debug...\"");
//
//

  typedef TrilinosWrappers::MPI::BlockVector VEC;

  Solver<VEC> solver(comm);
  DAETimeIntegrator<VEC> ode(solver);

  ParameterAcceptor::initialize(SOURCE_DIR "/parameters/sundials_interface_trilinos_01.prm", "ode.prm");

  shared_ptr<VEC> sol = solver.create_new_vector();
  sol->add(1);

  shared_ptr<VEC> sol_dot = solver.create_new_vector();

  IndexSet is = sol->locally_owned_elements();

  deallog << "[IS] N_dofs: " << is.size()
          << ", locally: " << is.n_elements() << std::endl;

  deallog << "Dofs: " << solver.n_dofs()
          << ", vec size: " << sol->size() << std::endl;

  unsigned int max_steps=10000;


  ode.start_ode(*sol, *sol_dot, max_steps);


  return (0);
}
