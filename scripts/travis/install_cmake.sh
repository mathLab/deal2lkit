#!/bin/sh

PRG=$PWD/programs
CASA=$PWD

# installing cmake
if [ ! -d $PRG/cmake ]
then
  echo "-------------------------------------->installing cmake"
  wget --no-check-certificate https://cmake.org/files/v3.5/cmake-3.5.0-rc3-Linux-x86_64.tar.gz
  tar xzf cmake-3.5.0-rc3-Linux-x86_64.tar.gz -C $PRG
  rm cmake-3.5.0-rc3-Linux-x86_64.tar.gz
  mv $PRG/cmake-* $PRG/cmake
fi
export PATH=$PRG/cmake/bin:$PATH


