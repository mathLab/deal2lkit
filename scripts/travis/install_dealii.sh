#!/bin/sh

PRG=$PWD/programs
CASA=$PWD


# dealii
if [ ! -d $PRG/dealii ]
then
  echo "-------------------------------------->installing dealii"
  DST_INST=$PRG/dealii 
  cd $PRG
  git clone https://github.com/dealii/dealii.git dealii-tmp
  cd dealii-tmp
  mkdir build
  cd build
  cmake \
    -G Ninja \
    -D CMAKE_INSTALL_PREFIX:PATH=$DST_INST \
    -D DEAL_II_WITH_MPI:BOOL=ON \
    -D DEAL_II_WITH_THREADS:BOOL=OFF \
    -D DEAL_II_WITH_P4EST:BOOL=ON \
    .. > $CASA/dealii_cmake.log 2>&1
  cat summary.log
  ninja -j3 -k0 ; ninja -j3 -k0; ninja -j2 -k0 ; ninja -j2 -k0; ninja -j1 -k0 ; ninja -j1
  ninja -j4 install > /dev/null
  cd $PRG
  rm -rf dealii-tmp
  cd $CASA
fi

