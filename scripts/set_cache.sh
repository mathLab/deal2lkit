#!/bin/sh

if [ ! -d programs ] 
then
	echo "create folder `$PWD`/programs"
	mkdir programs
fi

# installing cmake
if [ ! -d programs/cmake ]
then
	echo "installing cmake"
	wget http://www.cmake.org/files/v3.3/cmake-3.3.0-rc3-Linux-x86_64.tar.gz 
	tar xzf cmake-3.3.0-rc3-Linux-x86_64.tar.gz -C programs
	rm cmake-3.3.0-rc3-Linux-x86_64.tar.gz
	mv programs/cmake-* programs/cmake
fi

# installing ninja
if [ ! -d programs/ninja ]
then
	echo "installing ninja"
	cd programs
	git clone git://github.com/martine/ninja.git
	cd ninja
	git checkout release
	./configure.py --bootstrap
	cd ../..
fi

# dealii
if [ ! -d programs/dealii ]
then
	echo "installing dealii"
	mkdir programs/dealii
	export D2V=`git ls-remote --tags https://github.com/luca-heltai/dealii | tail -n 1 | awk -F"/" '{print $3}'`
	echo "Deal.II version $D2V"
	wget https://github.com/luca-heltai/dealii/releases/download/$D2V/dealii-travis-CI-build.tgz
	tar xfz dealii-travis-CI-build.tgz -C programs/dealii
fi

# astyle
if [ ! -d programs/astyle ]
then
	echo "Downloading and installing astyle."
	mkdir programs/astyle
	wget http://downloads.sourceforge.net/project/astyle/astyle/astyle%202.04/astyle_2.04_linux.tar.gz  > /dev/null
	tar xvfz astyle_2.04_linux.tar.gz -C programs
	cd programs/astyle/build/gcc
	make -j4 > /dev/null
	cd ../../../../
fi

