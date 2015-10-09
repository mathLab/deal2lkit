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

#ifndef _d2k_sak_data_h
#define _d2k_sak_data_h

#include <deal2lkit/config.h>
#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>

#include <boost/any.hpp>
#include <vector>
#include <algorithm>
#include <typeinfo>

#include <map>
#include <deal2lkit/utilities.h>

using namespace dealii;

/**
 * Store any amount of any type of data accessible by an identifier string.
 *
 * It is a std::map<std::string, boost::any>
 *
 * Maps are associative containers that store elements formed by
 * the combination of a key value and a mapped value, and which allows for
 * fast retrieval of individual elements based on their keys.
 *
 *
 * A typical usage of this class is the following
 * @code
 * SAKData data;
 * const unsigned int n_q = 5;
 * std::vector<double> v_double(n_q);
 * ...
 * data.add_copy<std::vector<double> >(v_double, "double_copy");
 * ...
 * auto &vd = data.get<std::vector<double> >("double_copy");

 * std::vector<int> v_int(n_q);
 * data.add_ref<std::vector<int> >(v_int, "int_ref");
 * v_int[0] = 7.0;
 * auto &vi = data.get<std::vector<int> >("int_ref");
 *
 * @endcode
 */

class SAKData : public Subscriptor
{
public:

  /**
   * @brief Add a copy of an object
   *
   * Add a copy of an object. The copied object is owned by the class.
   *
   */
  template <typename type>
  void add_copy (const type &entry, const std::string &name);

  /**
   * @brief Add a reference to an already existing object.
   *
   * Add a reference to an external object. The object is not owned by the class
   *
   */
  template <typename type>
  void add_ref (type &entry, const std::string &name);

  /**
   * @brief Access to stored data object by name.
   *
   * Find the object with given name, try to convert it to <tt>type</tt> and
   * return it. This function throws an exception if either the name does not
   * exist or if the conversion fails.
   *
   */
  template <typename type>
  type &get (const std::string &name);

  /**
   * @brief Read-only access to stored data object by name.
   *
   * Find the object with given name, try to convert it to <tt>type</tt> and
   * return it. This function throws an exception if either the name does not
   * exist or if the conversion fails.
   *
   */
  template <typename type>
  const type &get (const std::string &name) const;

  /**
   * @brief Query if we store the given object.
   *
   * Find out if we store an object with given name.
   *
   */
  inline bool have (const std::string &name) const;

  /**
   * @brief Print the name and type of the stored objects
   *
   * Print the name associated to each object and its type
   *
   */
  template <class STREAM>
  void print_info (STREAM &os);

  /// An entry with this name does not exist in the SAKData object.
  DeclException1(ExcNameNotFound, std::string,
                 << "No entry with the name " << arg1 << " exists.");

  /// The requested type and the stored type are different
  DeclException2(ExcTypeMismatch,
                 char *, char *,
                 << "The requested type " << arg1
                 << " and the stored type " << arg2
                 << " must coincide");


private:
  std::map<std::string, boost::any> mydata;
}; // end class

template <typename type>
void SAKData::add_copy (const type &entry, const std::string &name)
{
  mydata[name] = entry;
}

template <typename type>
void SAKData::add_ref (type &entry, const std::string &name)
{
  type *ptr = &entry;
  mydata[name] = ptr;
}

template <typename type>
type &SAKData::get(const std::string &name)
{
  Assert( mydata.find(name) != mydata.end(),
          ExcNameNotFound(name));

  type *p=NULL;

  if (mydata[name].type() == typeid(type *) )
    {
      p = boost::any_cast<type *>(mydata[name]);
    }
  else if (mydata[name].type() == typeid(type))
    {
      p = boost::any_cast<type>(&mydata[name]);
    }
  else
    {
      Assert(false,
             ExcTypeMismatch(typeid(type).name(),mydata[name].type().name()));
    }

  return *p;
}

template <typename type>
const type &SAKData::get(const std::string &name) const
{
  Assert( mydata.find(name) != mydata.end(),
          ExcNameNotFound(name));

  typedef std::map<std::string, boost::any>::const_iterator it_type;

  it_type it = mydata.find(name);

  if (it->second.type() == typeid(type *) )
    {
      const type *p = boost::any_cast<type *>(it->second);
      return *p;
    }
  else if (it->second.type() == typeid(type))
    {
      const type *p = boost::any_cast<type>(&it->second);
      return *p;
    }
  else
    {
      Assert(false,
             ExcTypeMismatch(typeid(type).name(),it->second.type().name()));
      const type *p=NULL;
      return *p;
    }
}


bool SAKData::have(const std::string &name) const
{
  return mydata.find(name) != mydata.end();
}

template <class STREAM>
inline
void SAKData::print_info(STREAM &os)
{
  for (std::map<std::string, boost::any>::iterator  it=mydata.begin(); it != mydata.end(); ++it)
    {
      os << it->first
         << '\t' << '\t'
         << demangle(it->second.type().name())
         << std::endl;
    }
}

#endif
