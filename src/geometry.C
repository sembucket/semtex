///////////////////////////////////////////////////////////////////////////////
// geometry.C: define geometrical properties for 2D quad X Fourier spaces.
//
// Most routines are inlined in header file Geometry.h
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <stdio.h>

#include <femdef.h>
#include <Utility.h>
#include <Geometry.h>
#include <Femlib.h>


integer            Geometry::np    = 0;
integer            Geometry::nz    = 0;
integer            Geometry::nel   = 0;
integer            Geometry::psize = 0;
Geometry::CoordSys Geometry::csys  = Geometry::Cartesian;


void Geometry::set (const integer  NP,
		    const integer  NZ,
		    const integer  NE,
		    const CoordSys CS)
// ---------------------------------------------------------------------------
// Load values of static internal variables.
//
// NB: the value of psize is the value of nPlane, but rounded up to the
// next even number in case nPlane is odd.  This is to simplify the handling
// of Fourier transforms, which can be based on a real--complex transform
// on some platforms.
// ---------------------------------------------------------------------------
{
  char routine[] = "Geometry::set", err[StrMax];

  np = NP; nz = NZ; nel = NE; csys = CS;

  psize = nPlane();
  if (psize & 1) psize++;

  if (nz > 1 && nz & 1) {
    sprintf (err, "nz must be even (%1d)", nz);
    message (routine, err, ERROR);
  }

#if (defined(VECFFT))
  int nfact, factors[32];
  if (nz > 1) {
    Femlib::primes23 (nz, nfact, factors);
    if (!nfact) {
      sprintf (err, "nz (%1d) must have prime factors 2 & 3", nz);
      message (routine, err, ERROR);
    }
#endif
}
