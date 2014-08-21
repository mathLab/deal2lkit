#ifndef UTILITIES_HLT
#define UTILITIES_HLT

#include <deal.II/base/utilities.h>
#include <deal.II/base/smartpointer.h>

using namespace dealii;

template <typename TYPE>
void smart_delete (SmartPointer<TYPE> &sp) {
  if(sp) {
    TYPE * p = sp;
    sp = 0;
    delete p;
  }
}

#endif
