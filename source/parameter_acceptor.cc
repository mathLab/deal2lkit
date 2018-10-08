//-----------------------------------------------------------
//
//    Copyright (C) 2015 - 2016 by the deal2lkit authors
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

#include <deal.II/base/path_search.h>
#include <deal.II/base/point.h>
#include <deal.II/base/revision.h>

#include <deal2lkit/parameter_acceptor.h>
#include <deal2lkit/revision.h>
#include <deal2lkit/utilities.h>

#include <fstream>

using namespace dealii;

D2K_NAMESPACE_OPEN

ParameterAcceptor::ParameterAcceptor(const std::string name)
  : dealii::ParameterAcceptor(name)
{}

D2K_NAMESPACE_CLOSE
