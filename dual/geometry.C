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

int_t              Geometry::pid   = 0;
int_t              Geometry::nproc = 0;
int_t              Geometry::ndim  = 0;
int_t              Geometry::np    = 0;
int_t              Geometry::nz    = 0;
int_t              Geometry::nzp   = 0;
int_t              Geometry::nel   = 0;
int_t              Geometry::psize = 0;
Geometry::CoordSys Geometry::csys  = Geometry::Cartesian;


void Geometry::set (const int_t  NP,
		    const int_t  NZ,
		    const int_t  NE,
		    const CoordSys CS)
// ---------------------------------------------------------------------------
// Load values of static internal variables.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set", err[StrMax];

  pid   = Femlib::ivalue ("I_PROC");
  nproc = Femlib::ivalue ("N_PROC");

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
