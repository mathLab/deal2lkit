#!/bin/sh

PRG=$PWD/programs
CASA=$PWD


# trilinos
if [ ! -d $PRG/trilinos -a ! -d $PRG/dealii ]
then
    $CASA/scripts/travis/install_trilinos.sh
else
    if [ ! -d $PRG/dealii ]
    then
	$CASA/scripts/travis/install_dealii.sh
    fi
fi

