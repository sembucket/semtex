///////////////////////////////////////////////////////////////////////////////
// filter.C: Legendre--Fourier lowpass filtering of field data.
//
// Copyright (c) 2000 Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "les.h"

static real       *FourierMask = 0, *LegendreMask = 0;
static real       *Du = 0, *Dt = 0;
static const real *Iu = 0, *It = 0;


void initFilters ()
// ---------------------------------------------------------------------------
// Set up 1D filters for smoothing data.  These are to be applied as
// spectral coefficient multipliers in fully-transformed space.  Also
// compute matrices required for Legendre forward transform+filter,
// inverse transform.
//
// Both sets of filters are of erfc (or Boyd--Vandeven) type.  The
// Fourier filter is of length NZ / 2 + 1 (so notionally it takes the
// Nyquist frequency into account too), and is zero-phase (i.e. only
// real coefficients).  The Legendre filter is long enough to fit
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
// where x = F for Fourier mask, L for Legendre mask.
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
  integer      i, n, nm, order;
  real         lag, atten;
  const real*  dlt;
  vector<real> work();

  if (FourierMask) return;	// -- Already initialised!

  // -- Fourier mask.

  n     = Geometry::nZ();
  nm    = Geometry::nMode();
  order = (integer) Femlib::value ("F_ORDER");
  lag   =           Femlib::value ("F_ROLL" );
  atten =           Femlib::value ("F_ATTEN");
  work.setSize (nm + 1);
  FourierMask = new real [n];

  Femlib::erfcFilter (nm, order, (integer) (lag * nm), atten, work());
  FourierMask[0] = work[ 0];
  Fouriermask[1] = work[nm];
  for (i = 1; i < nm; i++) FourierMask[2*i] = FourierMask[2*i+1] = work[i];

  // -- Legendre mask.
  
  n     = Geometry::nP();
  nm    = n - 1;
  order = (integer) Femlib::value ("L_ORDER");
  lag   =           Femlib::value ("L_ROLL" );
  atten =           Femlib::value ("L_ATTEN");

  LegendreMask = new real [n];
  Femlib::erfcFilter (nm, order, (integer) (lag * nm), atten, LegendreMask);

  // -- Legendre transform+filter (Du, Dt) and inversion (Iu, It) matrices.

  Du = new real [(size_t) 2 * n*n];
  Dt = Du + n*n;

  Femlib::legTran (n, &dlt, 0, &Iu, &It, 0, 0);
  for (i = 0; i < n; i++)
    Veclib::smul (n, LegendreMask[i], dlt + i*n, 1, Du + i*n, 1);
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      Dt[Veclib::row_major (j, i, n)] = Du[Veclib::row_major (i, j, n)];
}


void lowpass (real* data)
// ---------------------------------------------------------------------------
// Carry out lowpass filtering of data, which is assumed to be
// supplied in Fourier-transformed state, and to have usual AuxField
// structure.
// ---------------------------------------------------------------------------
{
  register integer i;
  const integer    pid = Geometry::ProcID();
  const integer    np  = Geometry::nP();
  const integer    nP  = Geometry::nPlane();
  const integer    nZP = Geometry::nZProc();
  vector<real>     tmp (nP);

  if (!FourierMask) initFilters();

  for (i = 0; i < nZP; i++) {
    Veclib::zero  (nP, tmp(), 1);
    Femlib::grad2 (np, data+i*nP, data+i*nP, tmp(), tmp(), Du, Dt, np, nel);
    Veclib::smul  (nP, FourierMask[i+pid*nZP], tmp(), 1, tmp(), 1);
    Veclib::zero  (nP, data+i*nP, 1);
    Femlib::grad2 (np, tmp(), tmp(), data+i*nP, data+i*nP, Iu, It, np, nel);
  }
}
