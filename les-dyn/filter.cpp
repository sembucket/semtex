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
// By default, both sets of filters are of erfc (or Boyd--Vandeven)
// type.  The Fourier filter is of length NZ / 2 + 1 (so notionally it
// takes the Nyquist frequency into account too), and is zero-phase
// (i.e. only real coefficients).  The polynomial filter is long
// enough to fit along the side of an element (i.e. Geometry::nP()
// long).
// 
// If PROJ is defined for compilation, the in-plane filtering takes
// place by a projection to a set of Lagrange interpolants of half the
// order of the basis selcted for the simulation.  In this case, no
// filter parameters are needed for the in-plane filtering.
//
// The Fourier mask is returned with length nZ, and has same structure
// as the planes of data (k=0.Re, k=Nyquist.Re, k=1.Re, k=1.Im,
// k=2.Re, ...).
//
// Each BVD filter is defined by three parameters (supplied in the
// TOKENS section)
//
//   x_ROLL,  real number in 0--1 that defines where filter roll-off begins;
//   x_ORDER, integer that defines filter order (e.g. 2 or more);
//   x_ATTEN, real number in 0--1 that defines degree of attenuation,
//
// where x = F for Fourier mask, P for polynomial mask.
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

  if (FourierMask) return;	// -- Already initialised.

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
    work.resize (nm + 1);

    Femlib::erfcFilter (nm, order, lag * nm, atten, &work[0]);
    FourierMask[0] = work[ 0];
    FourierMask[1] = work[nm];
    for (i = 1; i < nm; i++) FourierMask[2*i] = FourierMask[2*i+1] = work[i];
  }

  // -- polynomial mask.

#if defined (PROJ) // -- Projection to half-order Lagrange interpolants.

  const real **PF, **PT, **IB, **IT;

  n  = Geometry::nP();
  nm = (n + 1) / 2;
  Iu = new real [(size_t) (2 * n*n)];
  It = Iu + n*n;
  
  Femlib::mesh (GLL, GLL, n,  nm, 0, &PF, &PT, 0, 0);
  Femlib::mesh (GLL, GLL, nm, n,  0, &IB, &IT, 0, 0);
  Blas::mxm    (*IB, n, *PF, nm, Iu, n);
  Blas::mxm    (*PT, n, *IT, nm, It, n);

  ROOTONLY cout << "-- Spectral projection filter." << endl;

#else              // -- Parameterised Boyd--Vandeven filter.

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

#if 1
  ROOTONLY cout << "-- Modal DPT filter." << endl;
  Femlib::modTran (n, &Du, &Dt, 0, &dpt, 0, 0);
#else
  ROOTONLY cout << "-- Legendre DPT filter." << endl;
  Femlib::legTran (n, &Du, &Dt, 0, &dpt, 0, 0);
#endif

  for (i = 0; i < n; i++)
    Veclib::smul (n, PolyMask[i], dpt + i*n, 1, It + i*n, 1);

  Blas::mxm    (Dt, n, It, n, Iu, n);
  Veclib::copy (n*n,   Iu, 1, It, 1);

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      Iu[Veclib::row_major (j, i, n)] = It[Veclib::row_major (i, j, n)];

#endif
}


void lowpass (real* data)
// ---------------------------------------------------------------------------
// Carry out lowpass filtering of data, which is assumed to be
// supplied in Fourier-transformed state, and to have usual AuxField
// structure.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer pid = Geometry::procID();
  const integer np  = Geometry::nP();
  const integer nP  = Geometry::planeSize();
  const integer nZP = Geometry::nZProc();
  const integer nel = Geometry::nElmt();
  vector<real>  tmp (nP);
  real*         dataplane;

  if (!FourierMask) initFilters();

#if 1

  for (i = 0; i < nZP; i++) {
    dataplane = data + i*nP;
    Femlib::tpr2d (dataplane, dataplane, &tmp[0], Iu, It, np, np, nel);
    Veclib::smul  (nP, FourierMask[i+pid*nZP], dataplane, 1, dataplane, 1);
  }

#else

  for (i = 0; i < nZP; i++) {
    dataplane = data + i*nP;
    Femlib::tpr2d (dataplane, dataplane, &tmp[0], Du, Dt, np, np, nel);
    Femlib::tpr2d (dataplane, dataplane, &tmp[0], Iu, It, np, np, nel);
    Veclib::smul  (nP, FourierMask[i+pid*nZP], dataplane, 1, dataplane, 1);
  }

#endif
}













