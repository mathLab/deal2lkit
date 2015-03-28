#include "parameter_acceptor.h"
#include "utilities.h"


// Static empty class list
std::vector<SmartPointer<ParameterAcceptor> > ParameterAcceptor::class_list;

ParameterAcceptor::ParameterAcceptor(const std::string name) :
    acceptor_id(class_list.size()),
    section_name(name)
{
    SmartPointer<ParameterAcceptor> pt(this, type(*this).c_str());
    class_list.push_back(pt);
}


ParameterAcceptor::~ParameterAcceptor() {
    class_list[acceptor_id] = 0;
}

std::string ParameterAcceptor::get_section_name() const {
    return (section_name != "" ? section_name : type(*this));
}


void ParameterAcceptor::parse_all_parameters(ParameterHandler &prm) {
    for(unsigned int i=0; i< class_list.size(); ++i) {
        prm.enter_subsection(class_list[i]->get_section_name());
        class_list[i]->parse_parameters(prm);
        prm.leave_subsection();
    }
}

void ParameterAcceptor::declare_all_parameters(ParameterHandler &prm) {
    for(unsigned int i=0; i< class_list.size(); ++i) {
        prm.enter_subsection(class_list[i]->get_section_name());
        class_list[i]->declare_parameters(prm);
        prm.leave_subsection();
    }
}

void ParameterAcceptor::parse_parameters(ParameterHandler &prm) {
    for(auto it = parameters.begin(); it != parameters.end(); ++it) {
        if(it->second.type() == typeid(std::string *)) {
            *(boost::any_cast<std::string*>(it->second)) = prm.get(it->first);
        }
        if(it->second.type() == typeid(double *)) {
            *(boost::any_cast<double*>(it->second)) = prm.get_double(it->first);
        }
        if(it->second.type() == typeid(int *)) {
            *(boost::any_cast<int*>(it->second)) = prm.get_integer(it->first);
        }
        if(it->second.type() == typeid(unsigned int *)) {
            *(boost::any_cast<unsigned int*>(it->second)) = prm.get_integer(it->first);
        }
        if(it->second.type() == typeid(bool *)) {
            *(boost::any_cast<bool*>(it->second)) = prm.get_bool(it->first);
        }
    }
}


void ParameterAcceptor::parse_parameters_call_back() {};

