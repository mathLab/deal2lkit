#!/bin/sh

PRG=$PWD/programs
CASA=$PWD

# astyle
if [ ! -d $PRG/astyle ]
then
  echo "Downloading and installing astyle."
  mkdir $PRG/astyle
  wget http://downloads.sourceforge.net/project/astyle/astyle/astyle%202.04/astyle_2.04_linux.tar.gz  > /dev/null
  tar xfz astyle_2.04_linux.tar.gz -C $PRG > /dev/null
  cd $PRG/astyle/build/gcc
  make -j4 > /dev/null
  cd $CASA
fi
