///////////////////////////////////////////////////////////////////////////////
// geometry.C: define geometrical properties for 2D quad X Fourier spaces.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// Most routines are inlined in header file geometry.h.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <cstdio>
#include <iostream>

#include <cfemdef.h>
#include <utility.h>
#include <femlib.h>
#include <geometry.h>

int_t Geometry::_pid    = 0;
int_t Geometry::_nproc  = 0;
int_t Geometry::_np     = 0;
int_t Geometry::_nz     = 0;
int_t Geometry::_nzp    = 0;
int_t Geometry::_nel    = 0;
int_t Geometry::_psize  = 0;

int_t Geometry::_npert  = 0;
int_t Geometry::_nbase  = 0;
int_t Geometry::_nslice = 0;
Geometry::CoordSys Geometry::_csys = Geometry::Cartesian;
Geometry::Category Geometry::_cat  = Geometry::O2_3D_SYMM;

void Geometry::set (const int_t nel  ,
		    const int_t npert)
// ---------------------------------------------------------------------------
// Load values of static internal variables.  Session file should
// already have been dealt with.  As well as being specific to stability
// analysis, this version is written for serial execution.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set";
  char        err[StrMax];

  _pid    = Femlib::ivalue ("I_PROC");
  _nproc  = Femlib::ivalue ("N_PROC");

  _np     = Femlib::ivalue ("N_P");

  _nbase  = Femlib::ivalue ("N_BASE");
  _nslice = Femlib::ivalue ("N_SLICE");
  _csys   = Femlib::ivalue ("CYLINDRICAL") ? 
                     Geometry::Cylindrical : Geometry::Cartesian;
  _npert  = npert;
  _nel    = nel;
  _psize  = nPlane() + (nPlane() % 2);

  _nz = _nzp = Femlib::ivalue ("N_Z");

  if      (_nbase == 2 && _npert == 2 && _nz == 1) _cat = O2_2D;
  else if (_nbase == 2 && _npert == 3 && _nz == 1) _cat = O2_3D_SYMM;
  else if (_nbase == 2 && _npert == 3 && _nz == 2) _cat = O2_3D;
  else if (_nbase == 3 && _npert == 3 && _nz == 1) _cat = SO2_2D;
  else if (_nbase == 3 && _npert == 3 && _nz == 2) _cat = SO2_3D;
  else {
    sprintf (err, "illegal: N_BASE = %1d, N_PERT = %1d, N_Z = %1d",
	     _nbase, _npert, _nz); message (routine, err, ERROR);
  }

  // -- Other sanity checks.

  if (_nproc  > 1) message (routine, "serial execution only",          ERROR);
  if (_nslice < 1) message (routine, "N_SLICE must be set in session", ERROR);
}


const char* Geometry::symmetry ()
// ---------------------------------------------------------------------------
// Return a string value of the spatial symmetry category of problem.
// ---------------------------------------------------------------------------
{
  switch (_cat) {
  case O2_2D:      return "O(2), 2D"; break;
  case O2_3D_SYMM: return "O(2), 3D, standing wave"; break;
  case O2_3D:      return "O(2), 3D, travelling wave"; break;
  case SO2_2D:     return "SO(2), 2D"; break;
  case SO2_3D:     return "SO(2), 3D"; break;
  }
}

