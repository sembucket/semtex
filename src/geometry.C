///////////////////////////////////////////////////////////////////////////////
// geometry.C: define geometrical properties for 2D quad X Fourier spaces.
//
// Most routines are inlined in header file Geometry.h
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

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
//
// The number of processors is restricted: it must either be 1 or an
// even number, and it must be less than or equal to the number of
// planes / 2.  Furthermore, the number of planes on a processor must
// be even, unless NZ == 1.
//
// NB: the value of psize is the value of nPlane, but rounded up if
// necessary to be an even number and also an integer multiple of the
// number of processors.  The even number restriction is to simplify
// the handling of Fourier transforms, which can be based on a
// real--complex transform on some platforms.  The restriction to be
// an integer multiple of the number of processors is to simplify the
// structure of memory exchanges required for Fourier transforms.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set", err[StrMax];

  pid   = (integer) Femlib::value ("I_PROC");
  nproc = (integer) Femlib::value ("N_PROC");

  np   = NP; nz = NZ; nel = NE; csys = CS;
  nzp  = nz / nproc;
  ndim = (nz > 1) ? 3 : 2;

  if (nz > 1 && nz & 1) {	// -- 3D problems must have NZ even.
    sprintf (err, "N_Z must be even (%1d)", nz);
    message (routine, err, ERROR);
  }

  if (nproc > 1) {		// -- Concurrent execution restrictions.
    if (nproc & 1) {
      sprintf (err, "No. of processors must be even (%1d)",
	       nproc);
      message (routine, err, ERROR);
    }

    if (nproc << 1 > nz) {
      sprintf (err, "No. of processors (%1d) can at most be half N_Z (%1d)",
	       nproc, nz);
      message (routine, err, ERROR);
    }

    if (nz % (2 * nproc)) {
      sprintf (err, "No. of planes (%1d) per processor (%1d) must be even",
	       nz, nproc);
      message (routine, err, ERROR);
    }

    psize  = nPlane();
    psize += 2 * nproc - nPlane() % (2 * nproc);

  } else {

    psize = nPlane() + (nPlane() % 2);
  }
}
