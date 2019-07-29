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

#include <deal.II/base/parameter_acceptor.h>

#include <deal2lkit/config.h>


D2K_NAMESPACE_OPEN

/**
 * Moved to deal.II base class. Old interface is deprecated.
 */
class ParameterAcceptor : public dealii::ParameterAcceptor
{
public:
  ParameterAcceptor(const std::string &section_name = "")
    : dealii::ParameterAcceptor(section_name)
  {
    dealii::ParameterAcceptor::declare_parameters_call_back.connect(
      [&]() { this->declare_parameters_call_back(); });

    dealii::ParameterAcceptor::parse_parameters_call_back.connect(
      [&]() { this->parse_parameters_call_back(); });
  };

  /**
   * Add a parameter the given parameter list. A pointer to the
   * parameter is stored, so that every time the default
   * parse_parameters() function is called, this parameter is
   * updated with the value contained in prm.
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
    *parameter = dealii::Patterns::Tools::Convert<T>::to_value(default_value,
                                                               pattern.clone());
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
                const std::string &documentation = std::string(),
                const dealii::Patterns::PatternBase &pattern =
                  *dealii::Patterns::Tools::Convert<T>::to_pattern(),
                dealii::ParameterHandler &prm = dealii::ParameterAcceptor::prm)
  {
    dealii::ParameterAcceptor::add_parameter(
      entry, parameter, documentation, prm, pattern);
  }

  /**
   * Empty call back functions, compatible with deal2lkit implementation.
   */
  virtual void
  declare_parameters_call_back(){};

  /**
   * Empty call back functions, compatible with deal2lkit implementation.
   */
  virtual void
  parse_parameters_call_back(){};
};

D2K_NAMESPACE_CLOSE

#endif
