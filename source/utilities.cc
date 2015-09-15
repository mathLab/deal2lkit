#include <deal2lkit/utilities.h>
#include <vector>
#include <fstream>
#include <chrono>
#include <thread>

#ifdef DEAL_II_SAK_WITH_BOOST
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include "boost/filesystem.hpp"
#else
#include <stdlib.h>
#include <sys/stat.h>
#endif

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
#ifdef DEAL_II_SAK_WITH_BOOST
  return boost::filesystem::exists( file );
#else
  struct stat st;
  lstat(file.c_str(), &st);
  if (S_ISREG(st.st_mode))
    return true;
  else
    return false;
#endif
}

bool dir_exists(const std::string &dir)
{
#ifdef DEAL_II_SAK_WITH_BOOST
  return boost::filesystem::exists( dir );
#else
  struct stat st;
  lstat(dir.c_str(), &st);
  if (S_ISDIR(st.st_mode))
    return true;
  else
    return false;
#endif
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
  std::string name_cleaned = name;
  name_cleaned.erase(std::remove(name_cleaned.begin(),name_cleaned.end(),' '),name_cleaned.end());
#ifdef DEAL_II_SAK_WITH_BOOST
  boost::filesystem::create_directories("./" + name_cleaned + "/");
#else
  std::string cmd = "";
  cmd = "mkdir -p " + name_cleaned + "/ ;";
  std::system(cmd.c_str());
#endif
  return dir_exists(name_cleaned + "/");
}

bool copy_files(const std::string &files, const std::string &destination)
{
  if (files!="")
    {
      create_directory("./"+destination);
      bool result = true;
      std::vector<std::string> strs;
      std::string new_file;
      strs = dealii::Utilities::split_string_list(files, ' ');
      for (size_t i = 0; i < strs.size(); i++)
        {
          if (file_exists(strs[i]))
            {
              new_file = destination+"/"+strs[i];
#ifdef DEAL_II_SAK_WITH_BOOST
              boost::filesystem::copy_file(strs[i], new_file,
                                           boost::filesystem::copy_option::overwrite_if_exists);
#else
              std::string cmd = "cp " + strs[i] + " " +new_file + " ;";
              std::system( cmd.c_str() );
#endif
            }
          result &= file_exists(new_file);
        }
      return result;
    }
  else
    {
      return false;
    }
}

bool copy_file(const std::string &file, const std::string &new_file)
{
#ifdef DEAL_II_SAK_WITH_BOOST
  boost::filesystem::copy_file(file, new_file,
                               boost::filesystem::copy_option::overwrite_if_exists);
#else
  std::string cmd = "cp " + file + " " + new_file + " ;" ;
  std::system( cmd.c_str() );
#endif
  return file_exists(new_file);
}

bool rename_file(const std::string &file, const std::string &new_file)
{
#ifdef DEAL_II_SAK_WITH_BOOST
  boost::filesystem::rename(file, new_file);
#else
  std::string cmd = "mv " + file + " " + new_file +" ;";
  std::system( cmd.c_str() );
#endif
  return file_exists(new_file);
}
