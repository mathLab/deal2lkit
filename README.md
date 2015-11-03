# deal2lkit: A ToolKit library for deal.II (www.dealii.org) [![Build Status](https://travis-ci.org/mathLab/deal2lkit.svg)](https://travis-ci.org/mathLab/deal2lkit)

Copyright (C) 2014 -- 2015 by Luca Heltai (1) 

(1) Scuola Internazionale Superiore di Studi Avanzati
    E-mail: luca.heltai@sissa.it

Main Contributors:

https://github.com/mathlab/deal2lkit/graphs/contributors

This library contains a collection of useful C++ tools and classes
developed by the numerical analysis group at SISSA mathlab
(www.mathlab.sissa.it) in the last years, based on the `deal.II`
library (www.dealii.org) .  In this library we collect all the
developed tools which were not subject to licence restrictions. If you
find this collection useful, feel free to download, use it and suggest
pull requests!

The official distribution is on GitHub, and you can clone the
repository using

	git clone https://github.com/mathlab/deal2lkit.git

This file is subject to LGPL version 2.1 or later and may not be
distributed without copyright and license information. Please refer to
section 5 and 6 of this file for further information on this license.

A detailed documentation constructed using Doxygen is accessible here:

http://mathlab.github.io/deal2lkit

## 1. Deal.II Requirements: ##

Currently, the distributed source code requires the development
version of the deal.II library. It will not work with deal.II 8.3 and
previous versions, due to incompatibility issues.

2. Installation procedure:
==========================

Clone (or fork) the deal.II ToolKit using 

	git clone https://github.com/mathlab/deal2lkit.git

The library can then be compiled by running

	mkdir build
	cd build
	cmake -DDEAL_II_DIR=/path/to/deal.II ..
	make

and tested using 
    
	make test

You can modify the resulting `CMakeCache.txt` to set a different
installation path, or you can call cmake with the additional option
`-DCMAKE_INSTALL_PREFIX=/path/to/install/dir`.
	
3. Extensive documentation:
===========================

A browsable documentation of the library is available at

http://mathlab.github.io/deal2lkit

If the user has the program Doxygen installed, this documentation can
be compiled for offline browsing by building the `doc` target (only
available if Doxygen was found during configure time):

	make doc

In this case, the documentation will be accessible in the subdirectory

	doc/html

If the user wants the deal.II documentation to be inlined, then the
file http://www.dealii.org/developer/doxygen/deal.tag should be
downloaded to the building directory before making the documentation,
for example, using wget

	wget http://www.dealii.org/developer/doxygen/deal.tag
	make doc

5. Licence Informations
=======================

See the LICENCE file in this directory.
