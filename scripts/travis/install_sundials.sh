#!/bin/sh

PRG=$PWD/programs
CASA=$PWD



# sundials
if [ ! -d $PRG/sundials ]
then
  echo "-------------------------------------->installing sundials"
  DST_INST=$PRG/sundials
  cd $PRG
  wget http://pkgs.fedoraproject.org/repo/extras/sundials/sundials-2.6.2.tar.gz/3deeb0ede9f514184c6bd83ecab77d95/sundials-2.6.2.tar.gz
  tar xzf sundials-2.6.2.tar.gz
  rm sundials-2.6.2.tar.gz
  mv sundials-2.6.2 sundials-tmp
  cd sundials-tmp
  mkdir build
  cd build
#  export PATH=$PRG/cmake/bin:$PATH
  export PATH=$PWD/programs/ninja:$PATH
  cmake \
    -G Ninja \
    -D CMAKE_INSTALL_PREFIX:PATH=$DST_INST \
    -D CXX_ENABLE:BOOL=ON \
    -D MPI_ENABLE:BOOL=ON \
    -D MPI_CC:PATH=$PRG/mpich/bin/mpicc \
    -D MPI_CXX:PATH=$PRG/mpich/bin/mpicxx \
    ..
  ninja -j4 install
  cd $PRG
  rm -rf sundials-tmp
  cd $CASA
fi
export SUNDIALS_DIR=$PRG/sundials
