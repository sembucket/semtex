///////////////////////////////////////////////////////////////////////////////
// geometry.C: define geometrical properties for 2D quad X Fourier spaces.
//
// Copyright (C) 1994--2001 Hugh Blackburn.
//
// Most routines are inlined in header file Geometry.h.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream.h>

#include <femdef.h>
#include <Utility.h>
#include <Geometry.h>
#include <Femlib.h>

int Geometry::_pid    = 0;
int Geometry::_nproc  = 0;
int Geometry::_ndim   = 0;
int Geometry::_np     = 0;
int Geometry::_nz     = 0;
int Geometry::_nzp    = 0;
int Geometry::_nel    = 0;
int Geometry::_psize  = 0;
int Geometry::_kfund  = 0;
int Geometry::_npert  = 0;
int Geometry::_nbase  = 0;
int Geometry::_nslice = 0;
Geometry::CoordSys Geometry::_csys = Geometry::Cartesian;


void Geometry::set (const int nel  ,
		    const int npert)
// ---------------------------------------------------------------------------
// Load values of static internal variables.  Session file should
// already have been dealt with.  As well as being specific to stability
// analysis, this version is written for serial execution.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set";

  _pid       = static_cast<integer>(Femlib::value ("I_PROC"));
  _nproc     = static_cast<integer>(Femlib::value ("N_PROC"));
  _kfund     = static_cast<integer>(Femlib::value ("K_FUND"));
  _np        = static_cast<integer>(Femlib::value ("N_POLY"));

  _nbase     = static_cast<integer>(Femlib::value ("N_BASE"));
  _nslice    = static_cast<integer>(Femlib::value ("N_SLICE"));
  _csys      = (static_cast<integer>(Femlib::value ("CYLINDRICAL"))) ? 
                               Geometry::Cylindrical : Geometry::Cartesian;
  _npert     = npert;
  _nel       = nel;
  _psize     = nPlane();

#if 1
  _nz = _nzp = static_cast<integer>(Femlib::value ("N_Z"));
#else
  _nz = _nzp = (_nbase == 3 && _npert == 3) ? 2 : 1;
#endif

  _ndim = (_nbase == _npert && _nz == 1) ? 2 : 3;
  
  // -- Sanity checks.

  if (_nproc > 1)
    message (routine, "serial execution only",                          ERROR);
  if (_nbase < 2 || _nbase > 3)
    message (routine, "N_BASE must be set in session file",             ERROR);
  if (_nslice < 1)
    message (routine, "N_SLICE must be set in session file",            ERROR);
  if (_npert < 2 || _npert > 3)
    message (routine, "restart file has too many or few fields",        ERROR);
  if (static_cast<integer>(Femlib::value ("N_Z")) != _nz)
    message (routine, "declared value of N_Z clashes with requirement", ERROR);
}
