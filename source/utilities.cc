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
#include <vector>
#include <fstream>
#ifdef DEAL_II_WITH_CXX11
#include <thread>
#endif

#include <stdlib.h>
#include <sys/stat.h>


D2K_NAMESPACE_OPEN

namespace
{
  struct handle
  {
    char *p;
    handle(char *ptr) : p(ptr) { }
    ~handle()
    {
      delete p;
    }
  };
}

std::string demangle(const char *name)
{
  int status = -4; // some arbitrary value to eliminate the compiler warning
  handle result( abi::__cxa_demangle(name, NULL, NULL, &status) );
  return (status==0) ? result.p : name ;
}

#ifdef DEAL_II_WITH_CXX11
void TimeUtilities::sleep(unsigned int t)
{
  std::this_thread::sleep_for(std::chrono::milliseconds(t));
}

void TimeUtilities::get_start_time()
{
  AssertThrow(status == true,
              ExcMessage("Use get_end_time() before reuse get_start_time().") );

  t_start = std::chrono::high_resolution_clock::now();

  status = false;
}

void TimeUtilities::get_end_time()
{
  AssertThrow(status == false,
              ExcMessage("Use get_start_time() before get_end_time().") );

  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  times.push_back( time_span.count() );

  status = true;
}

int TimeUtilities::get_num_measures()
{
  return times.size();
}
#endif

void append_to_file(const std::string &in_file, const std::string &out_file)
{
  std::ifstream ifile(in_file,  std::ios::in);
  std::ofstream ofile(out_file, std::ios::out | std::ios::app);
  if (ifile.is_open())
    ofile << ifile.rdbuf();

  return;
}

bool file_exists(const std::string &file)
{
  struct stat st;
  return (stat (file.c_str(), &st) == 0);
}

bool dir_exists(const std::string &dir)
{
  struct stat st;
  return (stat (dir.c_str(), &st) == 0);
}

unsigned int get_next_available_index_directory_name(const std::string &base, int n_digits, unsigned int start, unsigned int index_max)
{
  if (start<index_max)
    {
      if ( dir_exists( base + dealii::Utilities::int_to_string (start, n_digits) ) )
        return get_next_available_index_directory_name(base, n_digits, ++start, index_max);
      else
        return start;
    }
  else
    return index_max;
}

std::string get_next_available_directory_name(const std::string &base, int n_digits, unsigned int start, unsigned int index_max)
{
  unsigned int index = get_next_available_index_directory_name(base, n_digits, start, index_max);
  return base + dealii::Utilities::int_to_string (index, n_digits);
}

bool create_directory(const std::string &name)
{
  Assert((std::find(name.begin(), name.end(), ' ') == name.end()),
         ExcMessage("Invalid name of directory."));
  std::string cmd = "mkdir -p " + name;
  std::system(cmd.c_str());
  return dir_exists(name);
}

bool copy_files(const std::string &files, const std::string &destination)
{
  create_directory("./"+destination);
  bool result = true;
  std::vector<std::string> strs;
  std::string new_file;
  strs = dealii::Utilities::split_string_list(files, ' ');
  for (size_t i = 0; i < strs.size(); i++)
    {
      Assert(file_exists(strs[i]), ExcMessage("Invalid name of file"));
      new_file = destination+"/"+strs[i];
      std::string cmd = "cp " + strs[i] + " " + new_file;
      std::system( cmd.c_str() );
      result &= file_exists(new_file);
    }
  return result;
}

bool copy_file(const std::string &file, const std::string &new_file)
{
  Assert(file_exists(file),ExcMessage("No such file or directory"));
  std::string cmd = "cp " + file + " " + new_file ;
  std::system( cmd.c_str() );
  return file_exists(new_file);
}

bool rename_file(const std::string &file, const std::string &new_file)
{
  Assert(file_exists(file),ExcMessage("No such file or directory"));
  std::string cmd = "mv " + file + " " + new_file;
  std::system( cmd.c_str() );
  return file_exists(new_file);
}


D2K_NAMESPACE_CLOSE

