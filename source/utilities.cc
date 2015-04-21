#include "../include/utilities.h"

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

