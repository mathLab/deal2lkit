// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal2lkit authors
//
// This file is part of the deal2lkit library.
//
// The deal2lkit library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

/**
@page changes_after_1_0_0 Changes after Version 1.0.0

<p>
This is the list of changes made after the release of \dk version
1.0.0. All entries are signed with the names of the authors.
</p>



<!-- ----------- INCOMPATIBILITIES ----------------- -->

<a name="incompatible"></a>
<h3 style="color:red">Incompatibilities</h3>

<p style="color:red">
Following are a few modifications to the library that unfortunately
are incompatible with previous versions of the library, but which we
deem necessary for the future maintainability of the
library. Unfortunately, some of these changes will require
modifications to application programs. We apologize for the
inconvenience this causes.
</p>

<ol> 

<li> Changed: The ParsedFunction, ParsedMappedFunctions, and
ParsedDirichletBCs used to have a template parameter specifying the
number of components. This has been removed and the number of
components is now specified at construction time. 
(Luca Heltai, 2015/11/26)
</li>

<li> Removed: The class ParsedFiniteElement has no more the system and
preconditioner couplings and related variables and methods.  
(Alberto Sartori, 2015/12/01) 
</li>


</ol>


<!-- ----------- GENERAL IMPROVEMENTS ----------------- -->

<a name="general"></a>
<h3>General</h3>

<ol>
  <li> 
  New: Add tests and clean for ParsedAMGPreconditioner and 
  ParsedJacobiPreconditioner.
  <br>
  (Mauro Bardelloni 2016/05/02)
  </li>
  <li> 
  New: The class ParsedILUPreconditioner has been implemented.
  <br>
  (Mauro Bardelloni 2016/05/02)
  </li>
  <li> 
  New: The class ParsedZeroAverageConstraints has been implemented.
  <br>
  (Alberto Sartori 2016/01/13)
  </li>


</ol>


<!-- ----------- SPECIFIC IMPROVEMENTS ----------------- -->

<a name="specific"></a>
<h3>Specific improvements</h3>


<ol>
  <li> 
  New: The function set_initial_time() has been implemented for
 IMEXStepper and IDAInterface.
  <br>
  (Alberto Sartori 2016/05/25)
  </li>


</ol>

*/
