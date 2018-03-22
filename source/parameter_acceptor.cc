#include <deal2lkit/parameter_acceptor.h>

D2K_NAMESPACE_OPEN

ParameterAcceptor::ParameterAcceptor(const std::string &section_name) :
  dealii::ParameterAcceptor(section_name)
{
  auto sig = [&]()
  {
    parse_parameters_call_back();
  };
  dealii::ParameterAcceptor::parse_parameters_call_back.connect(sig);
};



ParameterAcceptor::~ParameterAcceptor()
{}



void ParameterAcceptor::parse_parameters_call_back()
{}

D2K_NAMESPACE_CLOSE
