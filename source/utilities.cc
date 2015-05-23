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

std::string get_next_available_index_directory_name(const std::string &base, int n_digits)
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
return base + Utilities::int_to_string (index, n_digits);
}


bool create_directory(const std::string &name)
{
  std::string cmd = "";
  #ifdef DEAL_II_SAK_WITH_BOOST
          return boost::filesystem::create_directories(name);
  #else
          cmd = "mkdir -p " + name;
          if( int( std::system( cmd.c_str() ) == 0 ) )
          {
            return true;
          }
          else
          {
            return false;
          }
  #endif
}
