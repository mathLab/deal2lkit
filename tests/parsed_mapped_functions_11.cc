// test acts_on_this_id


#include "tests.h"
#include "utilities.h"
#include "parameter_acceptor.h"
#include "parsed_mapped_functions.h"


int main ()
{
  initlog();
  ParsedMappedFunctions<2,3> pmf("Mapped Functions", "u,u,p",
                                 "0=u % 1=1 % 6=u;p",
                                 "0=x;y;0 % 1=0;0;0 % 6=y*k;0;k","k=1");

  ParameterAcceptor::initialize();
  ParameterAcceptor::prm.log_parameters(deallog);

  unsigned int id_true = 0;
  unsigned int id_false = 7;

  deallog << "Mapped function acts on id 0 (expected true) --> "
          << pmf.acts_on_id(id_true)
          <<std::endl;
  deallog << "Mapped function acts on id 7 (expected false) --> "
          << pmf.acts_on_id(id_false)
          <<std::endl;

}
