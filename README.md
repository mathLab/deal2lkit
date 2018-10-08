# deal2lkit: A ToolKit library for deal.II (www.dealii.org) [![Build Status](https://travis-ci.org/mathLab/deal2lkit.svg)](https://travis-ci.org/mathLab/deal2lkit)

Copyright (C) 2014 -- 2018 by Luca Heltai (1) 

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

Please consider citing the following work:

    
    @article{SartoriGiulianiBardelloni-2015-a,
	Author = {Alberto Sartori and Nicola Giuliani and Mauro Bardelloni and Luca Heltai},
	Journal = {SoftwareX},
	Title = {{deal2lkit: a Toolkit Library for deal.II}},
	Note = {To appear},
	Year = {2018}}


The official distribution is on GitHub, and you can clone the
repository using

	git clone https://github.com/mathlab/deal2lkit.git

This file is subject to LGPL version 2.1 or later and may not be
distributed without copyright and license information. Please refer to
section 5 and 6 of this file for further information on this license.

A detailed documentation constructed using Doxygen is accessible here:

http://mathlab.github.io/deal2lkit

## Deal.II Requirements

deal2lkit versioning system follows closely deal.II one. Currently, 
deal.II 8.5.0, 9.0.0, and 9.1.0 are supported. Make sure you install
the same version of deal.II and deal2lkit, or you may get weird errors.

As new versions of deal.II are released, we provide a corresponding version
of deal2lkit that is compatible with the newly released deal.II version, 
with the same versioning number 


## Installation procedure:

Clone (or fork) the deal.II ToolKit using 

	git clone https://github.com/mathlab/deal2lkit.git

The library can then be compiled by running

	mkdir build
	cd build
	cmake -DDEAL_II_DIR=/path/to/deal.II ..
	make

and tested (if you enabled the cmake option `D2K_ENABLE_TESTING`) using 
    
	make setup_tests
	ctest

You can modify the resulting `CMakeCache.txt` to set a different
installation path, or you can call cmake with the additional option
`-DCMAKE_INSTALL_PREFIX=/path/to/install/dir`.
	
## Extensive documentation:

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

## Citing this work

If you use this library, please consider citing the following works:


	@article{SartoriGiulianiBardelloni-2015-a,
	Author = {Alberto Sartori and Nicola Giuliani and Mauro Bardelloni and Luca Heltai},
	Journal = {SoftwareX},
	Title = {{deal2lkit: a Toolkit Library for deal.II}},
	Year = {2018}}

	
	@article{MaierBardelloniHeltai-2016-a,
	Author = {Matthias Maier and Mauro Bardelloni and Luca Heltai},
	Journal = {Computers and Mathematics with Applications},
	Number = {1},
	Pages = {1-24},
	Title = {{\texttt{LinearOperator} -- a generic, high-level expression syntax for linear algebra}},
	Volume = {72},
	Year = {2016}}

	
	@article{AlzettaArndtBangerth-2018-a,
	Author = {Giovanni Alzetta and Daniel Arndt and Wolfgang Bangerth and Vishal Boddu and Benjamin Brands and Denis Davydov and Rene Gassm{\"o}ller and Timo Heister and Luca Heltai and Katharina Kormann and Martin Kronbichler and Matthias Maier and Jean-Paul Pelteret and Bruno Turcksin and David Wells},
	Journal = {Journal of Numerical Mathematics},
	Title = {The deal.II Library, Version 9.0},
	Year = {2018}}

## Licence Information

See the LICENCE file in this directory.
