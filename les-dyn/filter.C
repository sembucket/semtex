///////////////////////////////////////////////////////////////////////////////
// filter.C: Polynomial--Fourier lowpass filtering of field data.
//
// Copyright (c) 2000 Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "les.h"

static real       *FourierMask = 0, *PolyMask = 0;
static const real *Du = 0, *Dt = 0;
static real       *Iu = 0, *It = 0;


void initFilters ()
// ---------------------------------------------------------------------------
// Set up 1D filters for smoothing data.  These are to be applied as
// spectral coefficient multipliers in fully-transformed space.  Also
// compute matrices required for polynomial forward transform+filter,
// inverse transform.
//
// Both sets of filters are of erfc (or Boyd--Vandeven) type.  The
// Fourier filter is of length NZ / 2 + 1 (so notionally it takes the
// Nyquist frequency into account too), and is zero-phase (i.e. only
// real coefficients).  The polynomial filter is long enough to fit
// along the side of an element (i.e. Geometry::nP() long).
//
// The Fourier mask is returned with length nZ, and has same structure
// as the planes of data (k=0.Re, k=Nyquist.Re, k=1.Re, k=1.Im,
// k=2.Re, ...).
//
// Each filter is defined by three parameters (supplied in the TOKENS section)
//   x_ROLL,  real number in 0--1 that defines where filter roll-off begins;
//   x_ORDER, integer that defines filter order (e.g. 2 or more);
//   x_ATTEN, real number in 0--1 that defines degree of attenuation,
// where x = F for Fourier mask, L for polynomial mask.
//
// @Article{lih97,
//  author = 	 {Julia G. Levin and Mohammed Iskandarani and
//                  Dale B. Haidvogel},
//  title = 	 {A Spectral Filtering Procedure for Eddy-Resolving
//                  Simulations with a Spectral Element Ocean Model},
//  journal = 	 JCP,
//  year = 	 1997,
//  volume =	 137,
//  pages =	 {130--154}
// }
// ---------------------------------------------------------------------------
{
  static char routine[]= "initFilters";
  integer      i, j, n, nm, order;
  real         lag, atten;
  const real*  dpt;
  vector<real> work;

  if (FourierMask) return;	// -- Already initialised!

  // -- Fourier mask.

  n     = Geometry::nZ();
  nm    = Geometry::nMode();
  order = (integer) Femlib::value ("F_ORDER");
  lag   =           Femlib::value ("F_ROLL" );
  atten =           Femlib::value ("F_ATTEN");

  FourierMask = new real [n];

  if (order < 1) {
    message (routine, "Filter properties not defined", WARNING);
    FourierMask[0] = 1.0;
    FourierMask[1] = 0.0;
    for (i = 1; i < nm; i++) FourierMask[2*i] = FourierMask[2*i+1] = 1.0;
  } else {
    work.setSize (nm + 1);

    Femlib::erfcFilter (nm, order, lag * nm, atten, work());
    FourierMask[0] = work[ 0];
    FourierMask[1] = work[nm];
    for (i = 1; i < nm; i++) FourierMask[2*i] = FourierMask[2*i+1] = work[i];
  }

  // -- polynomial mask.
  
  n     = Geometry::nP();
  nm    = n - 1;
  order = (integer) Femlib::value ("P_ORDER");
  lag   =           Femlib::value ("P_ROLL" );
  atten =           Femlib::value ("P_ATTEN");

  PolyMask = new real [n];

  if (order < 1) {
    message (routine, "Filter properties not defined", WARNING);
    for (i = 0; i < n; i++) PolyMask[i] = 1.0;
  } else
    Femlib::erfcFilter (nm, order, lag, atten, PolyMask);

  // -- polynomial transform+filter (Du, Dt) and inversion (Iu, It) matrices.

  Iu = new real [(size_t) (2 * n*n)];
  It = Iu + n*n;

  Femlib::modTran (n, &Du, &Dt, 0, &dpt, 0, 0);
  for (i = 0; i < n; i++)
    Veclib::smul (n, PolyMask[i], dpt + i*n, 1, It + i*n, 1);
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      Iu[Veclib::row_major (j, i, n)] = It[Veclib::row_major (i, j, n)];
}


void lowpass (real* data)
// ---------------------------------------------------------------------------
// Carry out lowpass filtering of data, which is assumed to be
// supplied in Fourier-transformed state, and to have usual AuxField
// structure.
// ---------------------------------------------------------------------------
{
#if 1
  integer       i;
  const integer pid = Geometry::procID();
  const integer np  = Geometry::nP();
  const integer nP  = Geometry::planeSize();
  const integer nZP = Geometry::nZProc();
  const integer nel = Geometry::nElmt();
  vector<real>  tmp (nP);
  real*         dataplane;

  if (!FourierMask) initFilters();

  for (i = 0; i < nZP; i++) {
    dataplane = data + i*nP;
    Femlib::tpr2d (dataplane, dataplane, tmp(), Du, Dt, np, nel);
    Femlib::tpr2d (dataplane, dataplane, tmp(), Iu, It, np, nel);
    Veclib::smul  (nP, FourierMask[i+pid*nZP], dataplane, 1, dataplane, 1);
  }
#else
  const integer    nP  = Geometry::planeSize();
  const integer    np  = Geometry::nP();
  const integer    nel = Geometry::nElmt();
  const integer    np2 = sqr(np);
  const integer    nZP = Geometry::nZProc();
  const integer    pid = Geometry::procID();

  register integer i, j, ij, p, q, pq;
  register real    P, Q;
  integer          n, k;
  vector<real>     work (np2);
  real             *dataplane, *src, *tmp = work();

  if (!FourierMask) initFilters();

  for (k = 0; k < nZP; k++) {
    for (src = data + k * nP, n = 0; n < nel; n++, src += np2) {
      Veclib::zero (np2, tmp, 1);
      for (ij = 0, i = 0; i < np; i++) {
	for (j = 0; j < np; j++, ij++) {
	  for (pq = 0, p = 0; p < np; p++) {
	    P = Du[Veclib::row_major (i, p, np)];
	    for (q = 0; q < np; q++, pq++) {
	      Q = Dt[Veclib::row_major (q, j, np)];
	      tmp[ij] += P * Q * src[pq];
	    }
	  }
	}
      }
      Veclib::copy (np2, tmp, 1, src, 1);
    }
  }

  for (k = 0; k < nZP; k++) {
    dataplane = data + i*nP;
    Veclib::smul  (nP, FourierMask[i+pid*nZP], dataplane, 1, dataplane, 1);
  }

  for (k = 0; k < nZP; k++) {
    for (src = data + k * nP, n = 0; n < nel; n++, src += np2) {
      Veclib::zero (np2, tmp, 1);
      for (ij = 0, i = 0; i < np; i++) {
	for (j = 0; j < np; j++, ij++) {
	  for (pq = 0, p = 0; p < np; p++) {
	    P = Iu[Veclib::row_major (i, p, np)];
	    for (q = 0; q < np; q++, pq++) {
	      Q = It[Veclib::row_major (q, j, np)];
	      tmp[ij] += P * Q * src[pq];
	    }
	  }
	}
      }
      Veclib::copy (np2, tmp, 1, src, 1);
    }
  }

#endif
}
