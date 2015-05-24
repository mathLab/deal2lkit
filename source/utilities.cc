#include "../include/utilities.h"

#ifdef DEAL_II_SAK_WITH_BOOST
#include "boost/filesystem.hpp"
#else
#include <stdlib.h>
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

// Nothing to do here, yet...

int get_next_available_index_directory_name(const std::string &base, int n_digits)
{
  unsigned int index = 0;
  std::string cmd = "";
#ifdef DEAL_II_SAK_WITH_BOOST
  while ( boost::filesystem::exists( base + Utilities::int_to_string (index, n_digits)) ) index++;
#else
  cmd = "test -d " + base + Utilities::int_to_string (index, n_digits);
  while ( int(std::system( cmd.c_str() )) == 0 )
    {
      index++;
      cmd = "test -d " + base + Utilities::int_to_string (index, n_digits);
    }
#endif
  return index;
}

std::string get_next_available_directory_name(const std::string &base, int n_digits)
{
  unsigned int index = get_next_available_index_directory_name(base, n_digits);
  return base + Utilities::int_to_string (index, n_digits);
}

bool exist_file(const std::string &file)
{
#ifdef DEAL_II_SAK_WITH_BOOST
  return boost::filesystem::exists( file );
#else
  std::string cmd = "";
  cmd = "test -f " + file;
  if ( int(std::system( cmd.c_str() )) == 0 )
    {
      return true;
    }
  else
    {
      return false;
    }
#endif
}


bool create_directory(const std::string &name)
{
  std::string name_cleaned = name;
  name_cleaned.erase(std::remove(name_cleaned.begin(),name_cleaned.end(),' '),name_cleaned.end());
#ifdef DEAL_II_SAK_WITH_BOOST
  return boost::filesystem::create_directories(name_cleaned + "/");
#else
  std::string cmd = "";
  cmd = "mkdir -p " + name_cleaned + "/";
  if ( int( std::system( cmd.c_str() ) == 0 ) )
    {
      return true;
    }
  else
    {
      return false;
    }
#endif
}

bool copy_files(const std::string &files, const std::string &destination)
{
  bool result = true;
#ifdef DEAL_II_SAK_WITH_BOOST
  if (exists(files) && files!="")
    {
      vector<string> strs;
      boost::split(strs,files,boost::is_any_of(" "));
      for (size_t i = 0; i < strs.size(); i++)
        result &= copy_file(strs[i],destination+"/"+files,
                            copy_option::overwrite_if_exists);
    }
#else
  std::string cmd1 = "for f in " + files + "; do test -e $f ; done";
  std::string cmd2 = "for f in " + files + "; do cp $f " + destination + "; done";
  if (int(std::system( cmd1.c_str() )) == 0 && files!="")
    {
      if (int(std::system( cmd2.c_str() )))
        {
          result &= true;
        }
      else
        {
          result &= false;
        }
    }
#endif
return result;
}
