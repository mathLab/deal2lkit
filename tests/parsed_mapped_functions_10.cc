// testing default constructor
// if no expressions are given
// ZeroFunction<spacedim>(n_components)
// is set on the specified ids


#include "tests.h"
#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_mapped_functions.h>


int main ()
{
  initlog();
  ParsedMappedFunctions<3,1> pmf("Mapped functions","","5=ALL % 3=ALL");

  ParameterAcceptor::initialize(SOURCE_DIR "/parameters/parsed_mapped_functions_10.prm", "used_parameters.prm");
  ParameterAcceptor::prm.log_parameters(deallog);

  Point<3> p(2,3,5);

  std::vector<unsigned int> ids = pmf.get_mapped_ids();
  for (unsigned int i=0; i<ids.size(); ++i)
    {
      deallog << "Component mask id " << ids[i] <<" "<<  pmf.get_mapped_mask(ids[i]) << std::endl;
      deallog << "Parsed Function on id: "<< ids[i] <<" "<<  (*pmf.get_mapped_function(ids[i])).value(p) << std::endl;
    }
}
