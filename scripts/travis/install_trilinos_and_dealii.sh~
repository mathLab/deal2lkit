#!/bin/sh

PRG=$PWD/programs
CASA=$PWD


# trilinos
if [ ! -d $PRG/trilinos ]
then
  echo "-------------------------------------->installing trilinos"
  DST_INST=$PRG/trilinos
  export PATH=$PRG/cmake/bin:$PATH
  cd $PRG
  git clone https://github.com/trilinos/trilinos.git trilinos-tmp
  cd trilinos-tmp
  mkdir build
  cd build
  cmake \
    -G Ninja \
    -D CMAKE_BUILD_TYPE:STRING=RELEASE \
    -D TPL_ENABLE_Boost:BOOL=OFF \
    -D TPL_ENABLE_BoostLib:BOOL=OFF \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D TPL_ENABLE_Netcdf:BOOL=OFF \
    -D CMAKE_INSTALL_PREFIX:PATH=$DST_INST \
    -D Trilinos_ENABLE_OpenMP:BOOL=OFF \
    -D BUILD_SHARED_LIBS:BOOL=ON \
    -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=FALSE \
    -D Trilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF \
    -D Trilinos_ENABLE_TESTS:BOOL=OFF \
    -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
    -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
    -D Trilinos_ENABLE_Epetra:BOOL=ON \
    -D Trilinos_ENABLE_NOX:BOOL=ON \
    -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
    -D Trilinos_ENABLE_Tpetra:BOOL=ON \
    -D Trilinos_ENABLE_Kokkos:BOOL=ON \
    -D Trilinos_ENABLE_Sacado:BOOL=ON \
    -D Trilinos_ENABLE_Amesos:BOOL=ON \
    -D Trilinos_ENABLE_AztecOO:BOOL=ON \
    -D Trilinos_ENABLE_Ifpack:BOOL=ON \
    -D Trilinos_ENABLE_Rythmos:BOOL=ON \
    -D Trilinos_ENABLE_Piro:BOOL=ON \
    -D Trilinos_ENABLE_MOOCHO:BOOL=ON \
    -D Trilinos_ENABLE_ML:BOOL=ON \
    -D Trilinos_ENABLE_MueLu:BOOL=ON \
    -D Trilinos_ENABLE_Komplex:BOOL=ON \
    -D Trilinos_ENABLE_Thyra:BOOL=ON \
    -D Trilinos_ENABLE_TrilinosCouplings:BOOL=ON \
    -D Trilinos_ENABLE_Fortran:BOOL=OFF \
    -D CMAKE_CXX_COMPILER:PATH=$PRG/mpich/bin/mpicxx\
    -D CMAKE_C_COMPILER:PATH=$PRG/mpich/bin/mpicc \
    .. > $CASA/trilinos_cmake.log 2>&1

  ninja -j3 
  ninja -j4 install > /dev/null
  cd $PRG
  rm -rf trilinos-tmp
  cd $CASA
fi

export TRILINOS_DIR=$PRG/trilinos
