#include "../include/utilities.h"
#include <vector>
#include <fstream>

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

int get_next_available_index_directory_name(const std::string &base, int n_digits)
{
  unsigned int index = 0;
  while ( dir_exists( base + dealii::Utilities::int_to_string (index, n_digits) ) ) index++;
  return index;
}

std::string get_next_available_directory_name(const std::string &base, int n_digits)
{
  unsigned int index = get_next_available_index_directory_name(base, n_digits);
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

fixed_lines::fixed_lines(int n, std::ostream &stream_out)
  :
  n_max(n),
  curr_line(0),
  stream_out(stream_out)
{}

int
fixed_lines::get_current_line()
{
  return curr_line;
}

void
fixed_lines::goto_previous_line(int n_line, bool erase)
{
  if (curr_line>0)
    for (unsigned int i = 0; i<n_line; ++i )
      {
        // go to the previous line:
        std::cout << "\e[A";
        // erase the line line:
        if (erase)
          std::cout << "\e[K";

        curr_line--;

        if (curr_line == 0)
          break;
      }
}

void
fixed_lines::print_line(std::string &txt, bool erase)
{
  if (curr_line == n_max)
    goto_previous_line(n_max);

  std::string erase_line = "\e[K";
  if (erase)
    std::cout << erase_line;
  std::cout << txt << std::endl;
  curr_line++;
}
