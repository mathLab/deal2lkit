#!/bin/sh

PRG=$PWD/programs
CASA=$PWD


# installing ninja
if [ ! -d $PRG/ninja ]
then
  echo "-------------------------------------->installing ninja"
  cd $PRG
  git clone git://github.com/martine/ninja.git
  cd ninja
  git checkout release
  ./configure.py --bootstrap > /dev/null
  cd $CASA
fi
export PATH=$PWD/programs/ninja:$PATH


