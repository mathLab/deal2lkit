#!/bin/sh

PRG=$PWD/programs
CASA=$PWD


# mpich
if [ ! -d $PRG/mpich ]
then
  echo "-------------------------------------->installing mpich"
  mkdir $PRG/mpich-tmp
  cd $PRG/mpich-tmp
  wget --no-check-certificate http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
  tar xvzf mpich-3.2.tar.gz > /dev/null
  mv mpich-3.2 mpich
   cd mpich
   mkdir build && cd build
   ../configure CFLAGS="-w" --prefix=$PRG/mpich --disable-static > /dev/null
   make -j4 install > /dev/null
   cd $PRG
   rm -rf $PRG/mpich-tmp
fi
export PATH=$PRG/mpich/bin:$PATH
