#!/bin/sh

PRG=$PWD/programs
CASA=$PWD


# p4est

if [ ! -d $PRG/p4est ]
then
    echo "-------------------------------------->installing p4est"
    DST_INST=$PRG/p4est
    cd $PRG
    git clone https://github.com/cburstedde/p4est.git p4est-tmp
    cd p4est-tmp
    git submodule init; git submodule update
    ./bootstrap > /dev/null
    ./configure --enable-mpi --prefix=$DST_INST > /dev/null
    make -j4 > /dev/null
    make install
    cd $CASA
    rm -rf $PRG/p4est-tmp
fi
export P4EST_DIR=$PRG/p4est
    

