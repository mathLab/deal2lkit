//-----------------------------------------------------------
//
//    Copyright (C) 2014 - 2016 by the deal2lkit authors
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

#ifndef d2k_parameter_acceptor_h
#define d2k_parameter_acceptor_h

#include <deal2lkit/config.h>
#include <deal2lkit/utilities.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_acceptor.h>
#include <boost/any.hpp>
#include <typeinfo>

using namespace dealii;

D2K_NAMESPACE_OPEN

/**
 * A deal2lkit backward-compatibe wrapper to dealii::ParameterAcceptor.
 */
class ParameterAcceptor : public dealii::ParameterAcceptor
{
public:
  /**
   * The constructor adds derived classes to the list of acceptors. If
   * a section name is specified, then this is used to scope the
   * parameters in the given section, otherwise a pretty printed
   * version of the derived class is used.
   */
  ParameterAcceptor(const std::string &section_name="");

  /**
   * The destructor sets to zero the pointer relative to this index,
   * so that it is safe to destroy the mother class.
   */
  virtual ~ParameterAcceptor();

  /**
   * Parse parameter call back. This function is called at the end
   * of parse_parameters, to allow users to process their parameters
   * right after they have been parsed. The default implementation
   * is empty.
   *
   * You can use this function, for example, to create a quadrature
   * rule after you have read how many quadrature points you wanted
   * to use from the parameter file.
   */
  virtual void parse_parameters_call_back();

  /**
   * Add a parameter the given parameter list. A pointer to the
   * parameter is stored, so that every time the default
   * parse_parameters() function is called, this parameter is
   * updated with the value contained in prm.
   */
  template <class T>
  void add_parameter(ParameterHandler &prm, T *parameter,
                     const std::string &entry,
                     const std::string &default_value,
                     const Patterns::PatternBase &pattern =
                       *Patterns::Tools::Convert<T>::to_pattern(),
                     const std::string &documentation=std::string())
  {
    AssertThrow(std::is_const<T>::value == false,
                ExcMessage("You tried to add a parameter using a const "
                           "variable. This is not allowed, since these "
                           "variables will be filled later on when "
                           "parsing the parameter."));
    *parameter = Patterns::Tools::Convert<T>::to_value(default_value, pattern.clone());
    prm.add_parameter(entry, *parameter, documentation, pattern);
  }


  /**
   * Add a parameter to the global parameter handler ParameterAcceptor::prm.
   * A pointer to the parameter is stored, so that every time the default
   * parse_parameters() function is called, this parameter is
   * updated with the value contained in prm. The default value of the
   * parameter is taken by the current content of the parameter, and the
   * default Pattern is constructed using the utilities in Patterns::Tools.
   *
   * Notice that this function has a slightly different behaviour with respect
   * to the other add_parameter() method, since it assumes that the global
   * parameter is always in its root section, and therefore before calling
   * prm.add_entry() it will enter in the sections specified by this class,
   * and leave all entered sections after having declared the variable, leaving
   * the parameter in the same state as before, but having inserted the entry
   * in the nested sections returned by get_section_path().
   */
  template <class T>
  void add_parameter(T &parameter,
                     const std::string &entry,
                     const std::string &documentation=std::string(),
                     const Patterns::PatternBase &pattern=
                       *Patterns::Tools::Convert<T>::to_pattern(),
                     ParameterHandler &prm=dealii::ParameterAcceptor::prm)
  {
    AssertThrow(std::is_const<T>::value == false,
                ExcMessage("You tried to add a parameter using a const "
                           "variable. This is not allowed, since these "
                           "variables will be filled later on when "
                           "parsing the parameter."));

    auto secs = get_section_path();
    for (auto sec : secs)
      prm.enter_subsection(sec);
    prm.add_parameter(entry, parameter, documentation, pattern);
    for (auto sec : secs)
      prm.leave_subsection();
  }
};

D2K_NAMESPACE_CLOSE

#endif

