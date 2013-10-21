///////////////////////////////////////////////////////////////////////////////
// geometry.C: define geometrical properties for 2D quad X Fourier spaces.
//
// NB: This is a modified version for dual: it is restricted to
// single-process operation and N_Z = 3.  In addition since no Fourier
// transformation will be used, there is no need to round up the plane
// size to be an even number.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// Most routines are inlined in header file Geometry.h
///////////////////////////////////////////////////////////////////////////////

static char CVS[] = "$Id$";

#include <cstdio>
#include <iostream>

using namespace std;

#include <cfemdef.h>
#include <utility.h>
#include <geometry.h>
#include <femlib.h>

int_t              Geometry::_pid   = 0;
int_t              Geometry::_nproc = 0;
int_t              Geometry::_ndim  = 0;
int_t              Geometry::_np    = 0;
int_t              Geometry::_nz    = 0;
int_t              Geometry::_nzp   = 0;
int_t              Geometry::_nel   = 0;
int_t              Geometry::_psize = 0;
Geometry::CoordSys Geometry::_csys  = Geometry::Cartesian;


void Geometry::set (const int_t  NP,
		    const int_t  NZ,
		    const int_t  NE,
		    const CoordSys CS)
// ---------------------------------------------------------------------------
// Load values of static internal variables.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set", err[StrMax];

  _pid   = Femlib::ivalue ("I_PROC");
  _nproc = Femlib::ivalue ("N_PROC");

  _np   = NP; _nz = NZ; _nel = NE; _csys = CS;
  _nzp  = _nz;
  _ndim = 3;

  if (_nz != 3) {
    sprintf (err, "dual needs N_Z == 3 (%1d)", _nz);
    message (routine, err, ERROR);
  }

  if (_nproc != 1) {
    sprintf (err, "This is a serial code (1 process, not %1d)", _nproc);
    message (routine, err, ERROR);
  }

  _psize = nPlane();
}
