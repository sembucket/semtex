///////////////////////////////////////////////////////////////////////////////
// geometry.C: define geometrical properties for 2D quad X Fourier spaces.
//
// NB: This is a modified version for dual: it is restricted to
// single-process operation and N_Z = 3.  In addition since no Fourier
// transformation will be used, there is no need to round up the plane
// size to be an even number.
//
// Copyright (C) 1994, 2000 Hugh Blackburn
//
// Most routines are inlined in header file Geometry.h
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream.h>

#include <femdef.h>
#include <Utility.h>
#include <Geometry.h>
#include <Femlib.h>

integer            Geometry::pid   = 0;
integer            Geometry::nproc = 0;
integer            Geometry::ndim  = 0;
integer            Geometry::np    = 0;
integer            Geometry::nz    = 0;
integer            Geometry::nzp   = 0;
integer            Geometry::nel   = 0;
integer            Geometry::psize = 0;
Geometry::CoordSys Geometry::csys  = Geometry::Cartesian;


void Geometry::set (const integer  NP,
		    const integer  NZ,
		    const integer  NE,
		    const CoordSys CS)
// ---------------------------------------------------------------------------
// Load values of static internal variables.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set", err[StrMax];

  pid   = (integer) Femlib::value ("I_PROC");
  nproc = (integer) Femlib::value ("N_PROC");

  np   = NP; nz = NZ; nel = NE; csys = CS;
  nzp  = nz;
  ndim = 3;

  if (nz != 3) {
    sprintf (err, "dual needs N_Z == 3 (%1d)", nz);
    message (routine, err, ERROR);
  }

  if (nproc != 1) {
    sprintf (err, "This is a serial code (1 process, not %1d)", nproc);
    message (routine, err, ERROR);
  }

  psize = nPlane();
}
