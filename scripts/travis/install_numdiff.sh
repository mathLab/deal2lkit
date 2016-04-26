#!/bin/sh

PRG=$PWD/programs
CASA=$PWD



# installing numdiff
if [ ! -d $PRG/numdiff ]
then
  echo "-------------------------------------->installing numdiff"
  mkdir $PRG/numdiff-tmp
  cd $PRG/numdiff-tmp
  wget http://mirror.lihnidos.org/GNU/savannah/numdiff/numdiff-5.8.1.tar.gz
  tar xzf numdiff-5.8.1.tar.gz
  rm numdiff-5.8.1.tar.gz
  cd numdiff-5.8.1
  DST_INST=$PRG/numdiff
  ./configure --prefix=$DST_INST > /dev/null
  make -j4 install > /dev/null
  cd $CASA
  rm -rf $PRG/numdiff-tmp
fi


