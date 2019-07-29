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

#include <deal2lkit/utilities.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <fstream>
#include <thread>
#include <vector>

using namespace dealii;

D2K_NAMESPACE_OPEN

void
append_to_file(const std::string &in_file, const std::string &out_file)
{
  std::ifstream ifile(in_file, std::ios::in);
  std::ofstream ofile(out_file, std::ios::out | std::ios::app);
  if (ifile.is_open())
    ofile << ifile.rdbuf();

  return;
}

bool
file_exists(const std::string &file)
{
  struct stat st;
  return (stat(file.c_str(), &st) == 0);
}

bool
dir_exists(const std::string &dir)
{
  struct stat st;
  return (stat(dir.c_str(), &st) == 0);
}

unsigned int
get_next_available_index_directory_name(const std::string &base,
                                        int                n_digits,
                                        unsigned int       start,
                                        unsigned int       index_max)
{
  if (start < index_max)
    {
      if (dir_exists(base + dealii::Utilities::int_to_string(start, n_digits)))
        return get_next_available_index_directory_name(base,
                                                       n_digits,
                                                       ++start,
                                                       index_max);
      else
        return start;
    }
  else
    return index_max;
}

std::string
get_next_available_directory_name(const std::string &base,
                                  int                n_digits,
                                  unsigned int       start,
                                  unsigned int       index_max)
{
  unsigned int index =
    get_next_available_index_directory_name(base, n_digits, start, index_max);
  return base + dealii::Utilities::int_to_string(index, n_digits);
}

bool
create_directory(const std::string &name)
{
  Assert((std::find(name.begin(), name.end(), ' ') == name.end()),
         ExcMessage("Invalid name of directory."));
  std::string cmd = "mkdir -p " + name;
  int         status;
  status = std::system(cmd.c_str());
  AssertThrow(status == 0, ExcCannottExecuteCommand(cmd));
  return dir_exists(name);
}

bool
copy_files(const std::string &files, const std::string &destination)
{
  create_directory("./" + destination);
  bool                     result = true;
  int                      status;
  std::vector<std::string> strs;
  std::string              new_file;
  strs = dealii::Utilities::split_string_list(files, ' ');
  for (size_t i = 0; i < strs.size(); i++)
    {
      Assert(file_exists(strs[i]), ExcMessage("Invalid name of file"));
      new_file        = destination + "/" + strs[i];
      std::string cmd = "cp " + strs[i] + " " + new_file;
      status          = std::system(cmd.c_str());
      AssertThrow(status == 0, ExcCannottExecuteCommand(cmd));
      result &= file_exists(new_file);
    }
  return result;
}

bool
copy_file(const std::string &file, const std::string &new_file)
{
  Assert(file_exists(file), ExcMessage("No such file or directory"));
  std::string cmd = "cp " + file + " " + new_file;
  int         status;
  status = std::system(cmd.c_str());
  AssertThrow(status == 0, ExcCannottExecuteCommand(cmd));
  return file_exists(new_file);
}

bool
rename_file(const std::string &file, const std::string &new_file)
{
  Assert(file_exists(file), ExcMessage("No such file or directory"));
  std::string cmd = "mv " + file + " " + new_file;
  int         status;
  status = std::system(cmd.c_str());
  AssertThrow(status == 0, ExcCannottExecuteCommand(cmd));
  return file_exists(new_file);
}

D2K_NAMESPACE_CLOSE
