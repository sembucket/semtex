///////////////////////////////////////////////////////////////////////////////
// geometry.C: define geometrical properties for 2D quad X Fourier spaces.
//
// Most routines are inlined in header file Geometry.h
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <stdio.h>

#include <Utility.h>
#include <Geometry.h>

int                Geometry::np   = 0;
int                Geometry::nz   = 0;
int                Geometry::nel  = 0;
Geometry::CoordSys Geometry::csys = Geometry::Cartesian;


void Geometry::set (const int      NP,
		    const int      NZ,
		    const int      NE,
		    const CoordSys CS)
// ---------------------------------------------------------------------------
// Load values of static internal variables.
// ---------------------------------------------------------------------------
{
  char routine[] = "Geometry::set", err[StrMax];

  np = NP; nz = NZ; nel = NE; csys = CS;

  if (nz > 1 && nz & 1) {
    sprintf (err, "nz must be even (%1d)", nz);
    message (routine, err, ERROR);
  }
}
