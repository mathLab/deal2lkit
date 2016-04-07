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
#include <thread>

#include <stdlib.h>
#include <sys/stat.h>


D2K_NAMESPACE_OPEN

struct handle
{
  char *p;
  handle(char *ptr) : p(ptr) { }
  ~handle()
  {
    delete p;
  }
};

std::string demangle(const char *name)
{
  int status = -4; // some arbitrary value to eliminate the compiler warning
  handle result( abi::__cxa_demangle(name, NULL, NULL, &status) );
  return (status==0) ? result.p : name ;
}

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
  int status;
  status = std::system(cmd.c_str());
  AssertThrow(status == 0, ExcCannottExecuteCommand(cmd));
  return dir_exists(name);
}

bool copy_files(const std::string &files, const std::string &destination)
{
  create_directory("./"+destination);
  bool result = true;
  int status;
  std::vector<std::string> strs;
  std::string new_file;
  strs = dealii::Utilities::split_string_list(files, ' ');
  for (size_t i = 0; i < strs.size(); i++)
    {
      Assert(file_exists(strs[i]), ExcMessage("Invalid name of file"));
      new_file = destination+"/"+strs[i];
      std::string cmd = "cp " + strs[i] + " " + new_file;
      status = std::system( cmd.c_str() );
      AssertThrow(status == 0, ExcCannottExecuteCommand(cmd));
      result &= file_exists(new_file);
    }
  return result;
}

bool copy_file(const std::string &file, const std::string &new_file)
{
  Assert(file_exists(file),ExcMessage("No such file or directory"));
  std::string cmd = "cp " + file + " " + new_file ;
  int status;
  status = std::system(cmd.c_str());
  AssertThrow(status == 0, ExcCannottExecuteCommand(cmd));
  return file_exists(new_file);
}

bool rename_file(const std::string &file, const std::string &new_file)
{
  Assert(file_exists(file),ExcMessage("No such file or directory"));
  std::string cmd = "mv " + file + " " + new_file;
  int status;
  status = std::system(cmd.c_str());
  AssertThrow(status == 0, ExcCannottExecuteCommand(cmd));
  return file_exists(new_file);
}


#ifdef D2K_WITH_SUNDIALS

#ifdef DEAL_II_WITH_MPI

#ifdef DEAL_II_WITH_TRILINOS

void copy(TrilinosWrappers::MPI::Vector &dst, const N_Vector &src)
{
  IndexSet is = dst.locally_owned_elements();
  AssertDimension(is.n_elements(), NV_LOCLENGTH_P(src));
  for (unsigned int i=0; i<is.n_elements(); ++i)
    {
      dst[is.nth_index_in_set(i)] = NV_Ith_P(src, i);
    }
  dst.compress(VectorOperation::insert);
}

void copy(N_Vector &dst, const TrilinosWrappers::MPI::Vector &src)
{
  IndexSet is = src.locally_owned_elements();
  AssertDimension(is.n_elements(), NV_LOCLENGTH_P(dst));
  for (unsigned int i=0; i<is.n_elements(); ++i)
    {
      NV_Ith_P(dst, i) = src[is.nth_index_in_set(i)];
    }
}

void copy(TrilinosWrappers::MPI::BlockVector &dst, const N_Vector &src)
{
  IndexSet is = dst.locally_owned_elements();
  AssertDimension(is.n_elements(), NV_LOCLENGTH_P(src));
  for (unsigned int i=0; i<is.n_elements(); ++i)
    {
      dst[is.nth_index_in_set(i)] = NV_Ith_P(src, i);
    }
  dst.compress(VectorOperation::insert);
}

void copy(N_Vector &dst, const TrilinosWrappers::MPI::BlockVector &src)
{
  IndexSet is = src.locally_owned_elements();
  AssertDimension(is.n_elements(), NV_LOCLENGTH_P(dst));
  for (unsigned int i=0; i<is.n_elements(); ++i)
    {
      NV_Ith_P(dst, i) = src[is.nth_index_in_set(i)];
    }
}

#endif //DEAL_II_WITH_TRILINOS

#ifdef DEAL_II_WITH_PETSC

void copy(PETScWrappers::MPI::Vector &dst, const N_Vector &src)
{
  IndexSet is = dst.locally_owned_elements();
  AssertDimension(is.n_elements(), NV_LOCLENGTH_P(src));
  for (unsigned int i=0; i<is.n_elements(); ++i)
    {
      dst[is.nth_index_in_set(i)] = NV_Ith_P(src, i);
    }
  dst.compress(VectorOperation::insert);
}

void copy(N_Vector &dst, const PETScWrappers::MPI::Vector &src)
{
  IndexSet is = src.locally_owned_elements();
  AssertDimension(is.n_elements(), NV_LOCLENGTH_P(dst));
  for (unsigned int i=0; i<is.n_elements(); ++i)
    {
      NV_Ith_P(dst, i) = src[is.nth_index_in_set(i)];
    }
}

void copy(PETScWrappers::MPI::BlockVector &dst, const N_Vector &src)
{
  IndexSet is = dst.locally_owned_elements();
  AssertDimension(is.n_elements(), NV_LOCLENGTH_P(src));
  for (unsigned int i=0; i<is.n_elements(); ++i)
    {
      dst[is.nth_index_in_set(i)] = NV_Ith_P(src, i);
    }
  dst.compress(VectorOperation::insert);
}

void copy(N_Vector &dst, const PETScWrappers::MPI::BlockVector &src)
{
  IndexSet is = src.locally_owned_elements();
  AssertDimension(is.n_elements(), NV_LOCLENGTH_P(dst));
  for (unsigned int i=0; i<is.n_elements(); ++i)
    {
      NV_Ith_P(dst, i) = src[is.nth_index_in_set(i)];
    }
}

#endif //DEAL_II_WITH_PETSC

#endif //mpi

void copy(BlockVector<double> &dst, const N_Vector &src)
{
#ifdef DEAL_II_WITH_MPI
  AssertDimension((unsigned int)NV_LOCLENGTH_P(src), dst.size());
#else
  AssertDimension((unsigned int)NV_LENGTH_S(src), dst.size());
#endif
  for (unsigned int i=0; i<dst.size(); ++i)
    {
#ifdef DEAL_II_WITH_MPI
      dst[i] = NV_Ith_P(src, i);
#else
      dst[i] = NV_Ith_S(src, i);
#endif
    }
}

void copy(N_Vector &dst, const BlockVector<double> &src)
{
#ifdef DEAL_II_WITH_MPI
  AssertDimension((unsigned int)NV_LOCLENGTH_P(dst), src.size());
#else
  AssertDimension((unsigned int)NV_LENGTH_S(dst), src.size());
#endif
  for (unsigned int i=0; i<src.size(); ++i)
    {
#ifdef DEAL_II_WITH_MPI
      NV_Ith_P(dst, i) = src[i];
#else
      NV_Ith_S(dst, i) = src[i];
#endif
    }
}

#endif // sundials



D2K_NAMESPACE_CLOSE

