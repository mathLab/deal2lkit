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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/smartpointer.h>

#include <boost/any.hpp>

#include <deal2lkit/config.h>
#include <deal2lkit/utilities.h>

#include <typeinfo>

D2K_NAMESPACE_OPEN

/**
 * A parameter acceptor base class. This class is here for backward
 * compatibility with deal2lkit version < 9.0.0.
 *
 * The ParameterAcceptor was ported to deal.II with some differences.
 *
 * See dealii::ParameterAcceptor for a detailed documentation.
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
  ParameterAcceptor(const std::string section_name = "");

  /**
   * The destructor sets to zero the pointer relative to this index,
   * so that it is safe to destroy the mother class.
   */
  virtual ~ParameterAcceptor() = default;

  /**
   * Add a parameter the given parameter list. This is a compatibility
   * wrapper to the new version in deal.II, where the order of the
   * input parameters is different.
   */
  template <class T>
  void
  add_parameter(
    dealii::ParameterHandler &           prm,
    T *                                  parameter,
    const std::string &                  entry,
    const std::string &                  default_value,
    const dealii::Patterns::PatternBase &pattern = dealii::Patterns::Anything(),
    const std::string &                  documentation = std::string())
  {
    *parameter = dealii::Patterns::Tools::Convert<T>::to_value(default_value);
    prm.add_parameter(entry, *parameter, documentation, pattern);
  }


  /**
   * Add a parameter to the global parameter handler ParameterAcceptor::prm.
   * A pointer to the parameter is stored, so that every time the default
   * parse_parameters() function is called, this parameter is
   * updated with the value contained in prm. The default value of the
   * parameter is taken by the current content of the parameter, and the
   * default Pattern is constructed using the method to_pattern().
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
  void
  add_parameter(T &                parameter,
                const std::string &entry,
                const std::string &documentation             = std::string(),
                const dealii::Patterns::PatternBase &pattern = *to_pattern(T()),
                dealii::ParameterHandler &prm = ParameterAcceptor::prm)
  {
    AssertThrow(std::is_const<T>::value == false,
                dealii::ExcMessage("You tried to add a parameter using a const "
                                   "variable. This is not allowed, since these "
                                   "variables will be filled later on when "
                                   "parsing the parameter."));

    enter_my_subsection(prm);
    prm.add_parameter(entry, parameter, documentation, pattern);
    leave_my_subsection(prm);
  }


  /**
   * Given a class T, construct its default pattern to be used when declaring
   * parameters.
   */
  template <class T>
  static std::shared_ptr<dealii::Patterns::PatternBase>
  to_pattern(const T &)
  {
    return SP(dealii::Patterns::Tools::Convert<T>::to_pattern().release());
  }

  /**
   * Given a string, fill the value of the given parameter.
   */
  template <class T>
  static T
  to_type(const std::string &s)
  {
    return dealii::Patterns::Tools::Convert<T>::to_type(s);
  }

  /**
   * Given a parameter, return a string containing the given parameter.
   */
  template <class T>
  static std::string
  to_string(const T &t)
  {
    return dealii::Patterns::Tools::Convert<T>::to_string(t);
  }
};


D2K_NAMESPACE_CLOSE

#endif
