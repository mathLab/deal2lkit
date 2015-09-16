#!/bin/bash

for d in `ls -d */`
do
	echo processing $d
	cd $d
	for f in `ls -p | grep -v /`
	do
		sed -i 's/2014/2014/' $f
		#git checkout -- $f
	done
	cd ..
done
for f in `ls -p | grep -v /`
do
	sed -i 's/2014/2014/' $f
done
