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


#ifndef d2k_parsed_function_object_h
#define d2k_parsed_function_object_h

#include <deal.II/base/std_cxx14/utility.h>

#include <deal2lkit/config.h>

D2K_NAMESPACE_OPEN

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

template <typename F, typename Tuple, size_t... S>
decltype(auto)
apply_tuple_impl(F &&fn, Tuple &&t, std::index_sequence<S...>)
{
  return std::forward<F>(fn)(std::get<S>(std::forward<Tuple>(t))...);
}


template <typename F, typename Tuple>
decltype(auto)
apply_from_tuple(F &&fn, Tuple &&t)
{
  std::size_t constexpr tSize =
    std::tuple_size<typename std::remove_reference<Tuple>::type>::value;
  return apply_tuple_impl(std::forward<F>(fn),
                          std::forward<Tuple>(t),
                          std_cxx14::make_index_sequence<tSize>());
}


template <class T, typename Tuple>

D2K_NAMESPACE_CLOSE

#endif
