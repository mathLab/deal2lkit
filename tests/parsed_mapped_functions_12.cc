// test basic functionalities
// loop over the mapped ids


#include "tests.h"
#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_mapped_functions.h>


int main ()
{
  initlog();
  ParsedMappedFunctions<2,3> pmf("Mapped Functions", "u,u,p",
                                 "0=u % 1=1 % 6=u;p",
                                 "0=t;y;0 % 1=t;0;0 % 6=t;0;k","k=1");

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  Point<2> p(2,3);

  std::vector<unsigned int> ids = pmf.get_mapped_ids();
  for (unsigned int i=0; i<ids.size(); ++i)
    {
      deallog << "Component mask id " << ids[i] <<" "<<  pmf.get_mapped_mask(ids[i]) << std::endl;
      deallog << "Parsed Function on id: "<< ids[i] <<" "<<  (*pmf.get_mapped_function(ids[i])).value(p) << std::endl;
    }
  const double newtime = 77.;
  pmf.set_time(newtime);
  for (unsigned int i=0; i<ids.size(); ++i)
    {
      deallog << "Component mask id " << ids[i] <<" "<<  pmf.get_mapped_mask(ids[i]) << std::endl;
      deallog << "Parsed Function on id: "<< ids[i] <<" "<<  (*pmf.get_mapped_function(ids[i])).value(p) << std::endl;
    }
}
