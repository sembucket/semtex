///////////////////////////////////////////////////////////////////////////////
// geometry.C: define geometrical properties for 2D quad X Fourier spaces.
//
// Copyright (C) 1994,2003 Hugh Blackburn.
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
Geometry::Category Geometry::_cat  = Geometry::O2_3D_SYMM;


void Geometry::set (const int nel  ,
		    const int npert)
// ---------------------------------------------------------------------------
// Load values of static internal variables.  Session file should
// already have been dealt with.  As well as being specific to stability
// analysis, this version is written for serial execution.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set";
  char        err[StrMax];

  _pid       = static_cast<int>(Femlib::value ("I_PROC"));
  _nproc     = static_cast<int>(Femlib::value ("N_PROC"));
  _kfund     = static_cast<int>(Femlib::value ("K_FUND"));
  _np        = static_cast<int>(Femlib::value ("N_POLY"));

  _nbase     = static_cast<int>(Femlib::value ("N_BASE"));
  _nslice    = static_cast<int>(Femlib::value ("N_SLICE"));
  _csys      = (static_cast<int>(Femlib::value ("CYLINDRICAL"))) ? 
                               Geometry::Cylindrical : Geometry::Cartesian;
  _npert     = npert;
  _nel       = nel;
  _psize     = nPlane() + (nPlane() % 2);

  _nz = _nzp = static_cast<int>(Femlib::value ("N_Z"));

//  _ndim = (_nbase == _npert && _nz == 1) ? 2 : 3;

  if      (_nbase == 2 && _npert == 2 && _nz == 1) _cat = O2_2D;
  else if (_nbase == 2 && _npert == 3 && _nz == 1) _cat = O2_3D_SYMM;
  else if (_nbase == 2 && _npert == 3 && _nz == 2) _cat = O2_3D;
  else if (_nbase == 3 && _npert == 3 && _nz == 2) _cat = SO2_3D;
  else {
    sprintf (err, "illegal: N_BASE = %1d, N_PERT = %1d, N_Z = %1d",
	     _nbase, _npert, _nz); message (routine, err, ERROR);
  }

  // -- Other sanity checks.

  if (_nproc  > 1) message (routine, "serial execution only",          ERROR);
  if (_nslice < 1) message (routine, "N_SLICE must be set in session", ERROR);
}
