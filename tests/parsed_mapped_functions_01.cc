// test basic functionalities


#include "tests.h"
#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/parsed_mapped_functions.h>


int main ()
{
  initlog();
  ParsedMappedFunctions<2,3> pmf("Mapped Functions", "u,u,p", "0=0;1 % 1=2 % 6=0;1;2","0=x;y;0 % 1=0;0;0 % 6=y*k;beta*y;k","k=1,beta=2");

  ParameterAcceptor::initialize(SOURCE_DIR "/parameters/parsed_mapped_functions_01.prm", "used_parameters.prm");
  ParameterAcceptor::prm.log_parameters(deallog);

  Point<2> p(2,3);

  std::vector<unsigned int> ids(3);
  for (unsigned int i=0; i<ids.size(); ++i)
    ids[i] = i;

  unsigned int id = 3;

  if (std::find(ids.begin(),ids.end(), id) != ids.end())
    deallog << "c'e" <<std::endl;


  deallog << "Component mask id 0: " <<  pmf.get_mapped_mask(0) << std::endl;
  deallog << "Component mask id 1: " <<  pmf.get_mapped_mask(1) << std::endl;
  deallog << "Component mask id 6: " <<  pmf.get_mapped_mask(6) << std::endl;
  deallog << "Parsed Function on id 0: " <<  (*pmf.get_mapped_function(0)).value(p) << std::endl;
  deallog << "Parsed Function on id 1: " <<  pmf.get_mapped_function(1)->value(p) << std::endl;
  deallog << "Parsed Function on id 6: " <<  pmf.get_mapped_function(6)->value(p) << std::endl;
}
