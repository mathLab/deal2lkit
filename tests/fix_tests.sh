#!/bin/bash

for d in `ls -d */`
do
	echo processing $d
	cd $d
	for f in `ls -p | grep -v /`
	do
		#sed -i 's/include "tests.h"/include "..\/tests.h"/' $f
		#git checkout -- $f
		mv $f ..
	done
	cd ..
	rm -rf $d
done
