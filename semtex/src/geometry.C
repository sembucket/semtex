///////////////////////////////////////////////////////////////////////////////
// geometry.C: define geometrical properties for 2D quad X Fourier spaces.
//
// Copyright (C) 1994--2001 Hugh Blackburn
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

integer            Geometry::_pid   = 0;
integer            Geometry::_nproc = 0;
integer            Geometry::_ndim  = 0;
integer            Geometry::_np    = 0;
integer            Geometry::_nz    = 0;
integer            Geometry::_nzp   = 0;
integer            Geometry::_nel   = 0;
integer            Geometry::_psize = 0;
integer            Geometry::_kfund = 0;
Geometry::CoordSys Geometry::_csys  = Geometry::Cartesian;


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
//
// With the introduction of Geometry::kFund(), the first non-zero
// Fourier mode for 3D problems can have an assigned integer
// wavenumber K_FUND, in addition to its wavelength parameter BETA.
// For cylindrical problems where axis BCs have a wavenumber
// dependence, we can then obtain the correct set of BCs even if the
// first represented non-zero mode does not have wavenumber=1.  In all
// other cases (Cartesian, or cylindrical problems that do not include
// the axis in the domain) K_FUND is irrelevant.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set", err[StrMax];

  _pid   = (integer) Femlib::value ("I_PROC");
  _nproc = (integer) Femlib::value ("N_PROC");
  _kfund = (integer) Femlib::value ("K_FUND");

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
