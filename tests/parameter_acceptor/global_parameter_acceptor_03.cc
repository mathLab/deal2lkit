//-----------------------------------------------------------
//
//    Copyright (C) 2015 by the deal2lkit authors
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
#include <deal.II/base/path_search.h>
#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>


using namespace deal2lkit;

class MyClass : public deal2lkit::ParameterAcceptor
{
public:
  MyClass(const std::string &name = "My Class",
          const unsigned int &i=0):
    deal2lkit::ParameterAcceptor(name+std::to_string(c++)),
    i(i)
  {
    add_parameter(this->i, "An integer");
  };

private:
  unsigned int i;
  static unsigned int c;
};

unsigned int MyClass::c=0;

class MyMasterClass : public deal2lkit::ParameterAcceptor
{
public:
  MyMasterClass(const std::string &name = "MyMasterClass",
                const double &d=0):
    deal2lkit::ParameterAcceptor(name),
    d(d),
    mc("Slave", 1)
  {
    add_parameter(this->d, "A double");
  };

private:
  double d;
  MyClass mc;
};


int main ()
{
  initlog();

  MyMasterClass  f00("/");
  MyMasterClass  f0("Relative");
  MyMasterClass  f1("RelativeTraling/");
  MyMasterClass  f2("/Absolute");
  MyMasterClass  f3("/AbsoluteTrailing/");
  MyMasterClass  f4("/AbsoluteNested/Class");
  MyMasterClass  f5("/AbsoluteNested/ClassTrailing/");

  std::string test = "/AbsoluteTrailing/";
  auto l = Utilities::split_string_list(test, '/');
  deallog << l.size() << std::endl;
  for (auto s: l)
    deallog << s << std::endl;

  deal2lkit::ParameterAcceptor::log_info();
  deal2lkit::ParameterAcceptor::prm.log_parameters(deallog);
}
