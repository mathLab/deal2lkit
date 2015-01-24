#include "parameter_acceptor.h"

// Static empty class list
std::vector<SmartPointer<ParameterAcceptor> > ParameterAcceptor::class_list;

ParameterAcceptor::ParameterAcceptor(const std::string name) :
  acceptor_id(class_list.size()),
  section_name(name)
{
  SmartPointer<ParameterAcceptor> pt(this, typeid(*this).name());
  class_list.push_back(pt);
}
  

ParameterAcceptor::~ParameterAcceptor() {
  class_list[acceptor_id] = 0;
}

