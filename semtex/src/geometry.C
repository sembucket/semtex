///////////////////////////////////////////////////////////////////////////////
// geometry.C: define geometrical properties for 2D quad X Fourier spaces.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// Most routines are inlined in header file geometry
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <cstdio>
#include <iostream>

#include <cfemdef.h>
#include <utility.h>
#include <geometry.h>
#include <femlib.h>

int_t Geometry::_pid   = 0;
int_t Geometry::_nproc = 0;
int_t Geometry::_ndim  = 0;
int_t Geometry::_np    = 0;
int_t Geometry::_nz    = 0;
int_t Geometry::_nzp   = 0;
int_t Geometry::_nel   = 0;
int_t Geometry::_psize = 0;
int_t Geometry::_kfund = 0;
Geometry::CoordSys Geometry::_csys  = Geometry::Cartesian;


void Geometry::set (const int_t    NP,
		    const int_t    NZ,
		    const int_t    NE,
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
// necessary to be an even number and also an int_t multiple of the
// number of processors.  The even number restriction is to simplify
// the handling of Fourier transforms, which can be based on a
// real--complex transform on some platforms.  The restriction to be
// an int_t multiple of the number of processors is to simplify the
// structure of memory exchanges required for Fourier transforms.
//
// With the introduction of Geometry::kFund(), the first non-zero
// Fourier mode for 3D problems can have an assigned int_t
// wavenumber K_FUND, in addition to its wavelength parameter BETA.
// For cylindrical problems where axis BCs have a wavenumber
// dependence, we can then obtain the correct set of BCs even if the
// first represented non-zero mode does not have wavenumber=1.  In all
// other cases (Cartesian, or cylindrical problems that do not include
// the axis in the domain) K_FUND is irrelevant.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set", err[StrMax];

  _pid   = Femlib::ivalue ("I_PROC");
  _nproc = Femlib::ivalue ("N_PROC");
  _kfund = Femlib::ivalue ("K_FUND");

  _np   = NP; _nz = NZ; _nel = NE; _csys = CS;
  _nzp  = _nz / _nproc;
  _ndim = (_nz > 1) ? 3 : 2;

  if (_nz > 1 && _nz & 1) {	// -- 3D problems must have NZ even.
    sprintf (err, "N_Z must be even (%1d)", _nz);
    message (routine, err, ERROR);
  }

  if (_nproc > 1) {		// -- Concurrent execution restrictions.
    if (_nproc & 1) {
      sprintf (err, "No. of processors must be even (%1d)",
	       _nproc);
      message (routine, err, ERROR);
    }

    if (_nproc << 1 > _nz) {
      sprintf (err, "No. of processors (%1d) can at most be half N_Z (%1d)",
	       _nproc, _nz);
      message (routine, err, ERROR);
    }

    if (_nz % (2 * _nproc)) {
      sprintf (err, "No. of planes (%1d) per processor (%1d) must be even",
	       _nz, _nproc);
      message (routine, err, ERROR);
    }

    _psize  = nPlane();
    _psize += 2 * _nproc - nPlane() % (2 * _nproc);

  } else {

    _psize = nPlane() + (nPlane() % 2);
  }
}
