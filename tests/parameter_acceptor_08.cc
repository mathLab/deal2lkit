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



#include "tests.h"
#include <deal.II/base/path_search.h>
#include <deal2lkit/utilities.h>
#include <deal2lkit/parameter_acceptor.h>


using namespace deal2lkit;

class FirstClass : public ParameterAcceptor
{
public:
  FirstClass(const std::string &name = "First Class"):
    ParameterAcceptor(name)
  {};

  virtual void declare_parameters(ParameterHandler &prm)
  {
    add_parameter(prm, &f_i, "First int", "3", Patterns::Integer(0));
    add_parameter(prm, &f_d, "First double", "7.7", Patterns::Double(0));
    add_parameter(prm, &f_b, "First bool", "true", Patterns::Bool());
    add_parameter(prm, &f_s, "First string", "hello");
  };

private:
  int         f_i;
  double      f_d;
  bool        f_b;
  std::string f_s;
};

class SecondClass : public ParameterAcceptor
{
public:
  SecondClass(const std::string &name = "Second Class"):
    ParameterAcceptor(name)
  {};

  virtual void declare_parameters(ParameterHandler &prm)
  {
    add_parameter(prm, &s_i, "Second int", "5", Patterns::Integer(0));
    add_parameter(prm, &s_d, "Second double", "9.9", Patterns::Double(0));
    add_parameter(prm, &s_b, "Second bool", "false", Patterns::Bool());
    add_parameter(prm, &s_s, "Second string", "bye bye");
  };

private:
  int         s_i;
  double      s_d;
  bool        s_b;
  std::string s_s;
};

int main ()
{
  initlog();

  FirstClass  f;
  SecondClass s;
  std::string output_name = "used_parameter_acceptor_08.prm";
  ParameterAcceptor::initialize(SOURCE_DIR "/parameters/parameter_acceptor_08.prm", output_name);
  ParameterAcceptor::prm.log_parameters(deallog);

}
