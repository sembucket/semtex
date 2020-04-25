# utility README {#utility_readme}

utility README
==============

Synopsis
--------
Utility programs for pre- and post-processing, coded in both C and C++.


Notes
-----

* See the brief descriptions provided in the utility Directory Reference.

* In general, the C programs provide actions which do not rely on
  spectral element shape function operators, but may e.g. do
  manipulations based on 2D elemental storage. Usually these are
  stand-alone pieces of code and do not link other object files to
  create an exceutable.

* On the other hand, the C++ programs can do operations which require
  knowledge of elemental shape functions etc, and which may build upon
  routines defined in src/*.cpp: quite often, the C++ programs also
  link object files compiled from src/*.cpp, and/or femlib, routines.

------------------------------------------------------------------------------

$Id: README.md,v 9.1 2020/01/06 04:35:44 hmb Exp $
