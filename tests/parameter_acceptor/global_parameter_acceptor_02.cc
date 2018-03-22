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

template<class T>
class TemplateClass : public deal2lkit::ParameterAcceptor
{
public:
  TemplateClass(const T &i, const std::string &var):
    deal2lkit::ParameterAcceptor("Class<"+var+" >"),
    t(i)
  {
    add_parameter(t, var);
  };

private:
  T t;
};

int main ()
{
  initlog();

  TemplateClass<std::string  > f00("ciao"         ,"std::string");
  TemplateClass<int          > f01(-1             ,"int");
  TemplateClass<unsigned int > f02(1              ,"unsigned int");
  TemplateClass<double       > f03(1e-4           ,"double");
  TemplateClass<bool         > f04(false          ,"bool");
  TemplateClass<Point<1>     > f05(Point<1>(0)    ,"Point<1>");
  TemplateClass<Point<2>     > f06(Point<2>(0,1)  ,"Point<2>");
  TemplateClass<Point<3>     > f07(Point<3>(0,1,2),"Point<3>");

  TemplateClass<std::vector<std::string  > > f10(std::vector<std::string  >(3, "ciao"          ),"std::vector<std::string>");
  TemplateClass<std::vector<int          > > f11(std::vector<int          >(3, -1              ),"std::vector<int>");
  TemplateClass<std::vector<unsigned int > > f12(std::vector<unsigned int >(3, 1               ),"std::vector<unsigned int>");
  TemplateClass<std::vector<double       > > f13(std::vector<double       >(3, 1e-4            ),"std::vector<double>");
  TemplateClass<std::vector<Point<1>     > > f15(std::vector<Point<1>     >(3, Point<1>(0)     ),"std::vector<Point<1>>");
  TemplateClass<std::vector<Point<2>     > > f16(std::vector<Point<2>     >(3, Point<2>(0,1)   ),"std::vector<Point<2>>");
  TemplateClass<std::vector<Point<3>     > > f17(std::vector<Point<3>     >(3, Point<3>(0,1,2) ),"std::vector<Point<3>>");

  TemplateClass<std::vector<std::vector<std::string  > > > f110(std::vector<std::vector<std::string  > >(2, std::vector<std::string  >(3, "ciao" )), "std::vector<std::vector<std::string>>");
  TemplateClass<std::vector<std::vector<int          > > > f111(std::vector<std::vector<int          > >(2, std::vector<int          >(3, -1     )), "std::vector<std::vector<int>>");
  TemplateClass<std::vector<std::vector<unsigned int > > > f112(std::vector<std::vector<unsigned int > >(2, std::vector<unsigned int >(3, 1      )), "std::vector<std::vector<unsigned int>>");
  TemplateClass<std::vector<std::vector<double       > > > f113(std::vector<std::vector<double       > >(2, std::vector<double       >(3, 1e-4   )), "std::vector<std::vector<double>>");

  std::system("rm -f parameters.prm");
  // The first round is here to create the file parameters.prm with default values
  deal2lkit::ParameterAcceptor::initialize("parameters.prm");
  // This is here to read the default values just created.
  deal2lkit::ParameterAcceptor::initialize("parameters.prm");
  deal2lkit::ParameterAcceptor::prm.log_parameters(deallog);
}
