/*****************************************************************************
 * energy.c: compute solution energy, enstrophy moments, other scalar
 * properties.
 *
 * Copyright (C) 1992, 1999 Hugh Blackburn
 * 
 * $Id$
 *****************************************************************************/

#include "iso.h"

#define EN(k1,k2,k3) ( MAG(U[1][(k1)][(k2)][(k3)]) + \
		       MAG(U[2][(k1)][(k2)][(k3)]) + \
		       MAG(U[3][(k1)][(k2)][(k3)]) )


real energyP (CVF            V   ,
	      const complex* Wtab)
/* ------------------------------------------------------------------------- *
 * Compute & return k = <UiUi>/2: diagnostic.  Do sums in PHYSICAL space.
 * Also: compare with energyF() below to check discrete Parseval identity.
 * ------------------------------------------------------------------------- */
{
  const int     ntot = N * N * N;
  register int  i, c;
  register real k;
  register real *u = &V[1][0][0][0].Re,
                *v = &V[2][0][0][0].Re,
                *w = &V[3][0][0][0].Re; 
  
  for (c = 1; c <= 3; c++)            /* --> Physical space. */
    rc3DFT (V[c], Wtab, INVERSE);

  k = 0.0;
  for (i = 0; i < ntot; i++) k += u[i]*u[i] + v[i]*v[i] + w[i]*w[i];
  k /= 2.0 * ntot;

  for (c = 1; c <= 3; c++) {          /* --> Fourier space. */
    rc3DFT  (V[c], Wtab, FORWARD);
    scaleFT (V[c]);
  }
  
  return k;
}


real energyF (const CVF U)
/* ------------------------------------------------------------------------- *
 * Compute & return k = <Ui Ui>/2.  Do sums in FOURIER space, using
 * Parseval's identity.
 * ------------------------------------------------------------------------- */
{
  const int        ntot = N * N * K;
  register int     c, i;
  register real    k;
  register complex *u;

  for (k = 0.0, c = 1; c <= 3; c++) {
    u = &U[c][0][0][0];
    for (i = 0; i < ntot; i++)
      k += MAG (u[i]);
  }

  return k;
}


real enstrophyF (const CVF U   ,
		 CVF       vort,
		 CF        work)
/* ------------------------------------------------------------------------- *
 * Compute & return k = <Omega_i Omega_i>/2.  Do sums in FOURIER space.
 * ------------------------------------------------------------------------- */
{
  const int        ntot = N * N * K;
  register int     c, i;
  register real    k;
  register complex *wi;

  curl (U, vort, work);

  for (k = 0.0, c = 1; c <= 3; c++) {
    wi = &vort[c][0][0][0];
    for (i = 0; i < ntot; i++)
      k += MAG(wi[i]);
  }

  return k;
}


real rmsEns (const CVF U)
/* ------------------------------------------------------------------------- *
 * Compute & return Omega, generalized enstrophy, of order 1.  This is just
 * the mean-squared derivative of the velocity field.  See Ref [5], eq. (3.1).
 * ------------------------------------------------------------------------- */
{
  register int  c, k1, b1, k2, b2, k3;
  register real omega;
  register real kSqrd;

  omega = 0.0;

  for (c = 1; c <= 3; c++)
    for (k1 = 1; k1 < K; k1++) {
      b1     = N - k1;
      kSqrd  = SQR (k1);
      omega += kSqrd * (MAG (U[c][k1][ 0][ 0]) +
			MAG (U[c][ 0][k1][ 0]) +
			MAG (U[c][ 0][ 0][k1]) );

      for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	b2     = N - k2;
	kSqrd  = SQR (k1) + SQR (k2);
	omega += kSqrd * (MAG (U[c][ 0][k1][k2]) + MAG (U[c][ 0][b1][k2]) +
			  MAG (U[c][k1][ 0][k2]) + MAG (U[c][b1][ 0][k2]) +
			  MAG (U[c][k1][k2][ 0]) + MAG (U[c][b1][k2][ 0]) );

	for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	  kSqrd  = SQR (k1) + SQR (k2) + SQR (k3);
	  omega += kSqrd * (MAG (U[c][k1][k2][k3]) + MAG (U[c][b1][k2][k3]) +
			    MAG (U[c][k1][b2][k3]) + MAG (U[c][b1][b2][k3]) );
	}
      }
    }
  
  return  omega;
}


real L2norm (const CF U)
/* ------------------------------------------------------------------------- *
 * Compute & return L2 norm = Sum U^2.  U supplied in FOURIER space.
 * ------------------------------------------------------------------------- */
{
  register int      i;
  register real     l2;
  register complex* u    = &U[0][0][0];
  const int         ntot = N * N * K;

  l2 = 0.0;
  for (i = 0; i < ntot; i++) l2 += MAG (u[i]);
  l2 *= 2.0;

  return l2;
}


real amaxF (const CF U)
/* ------------------------------------------------------------------------- *
 * Find the maximum value of scalar field U, given in PHYSICAL space.
 * ------------------------------------------------------------------------- */
{
  register int   i;
  register real  mx   = 0.0;
  register real* u    = &U[0][0][0].Re;
  const int      ntot = N * N * N;

  for (i = 0; i < ntot; i++) mx = MAX (fabs(u[i]), mx);

  return mx;
}


void normaliseVF (CVF IC)
/* ------------------------------------------------------------------------- *
 * Normalise velocity components to give k = 1.0.
 * 
 * IC components are supplied in FOURIER space.
 * ------------------------------------------------------------------------- */
{
  int        c;
  const real k = energyF (IC);
  
  for (c = 1; c <= 3; c++) scaleF (IC[c], 1.0 / sqrt (k));
}


void energySpec (const CVF U   ,
		 real*     spec)
/* ------------------------------------------------------------------------- *
 * Return the spherically-integrated energy spectrum E(k) in spec, which
 * is a vector K long.   Note that there will possibly be some energy
 * contributions at (rounded) wavenumbers that do not lie within K.
 *
 * Normalise spectrum so that its integral is q = kinetic energy.  Ignore
 * contributions on Nyquist planes.  
 * ------------------------------------------------------------------------- */
{
  const real   norm = 1.0;
  real         de;
  register int k, k1, b1, k2, b2, k3;

  memset (spec, '\0', K * sizeof (real));

  spec[0] += norm * EN (0, 0, 0);

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    de = norm * (EN (k1, 0, 0) + EN (0, k1, 0) + EN (0, 0, k1));
    spec[k1] += de;
    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;
      k  = (int) sqrt (k1*k1 + k2*k2);
      if (k >= K) continue;
      de = norm * (EN (0, k1, k2) + EN (0, b1, k2)); spec[k] += de;
      de = norm * (EN (k1, 0, k2) + EN (b1, 0, k2)); spec[k] += de;
      de = norm * (EN (k1, k2, 0) + EN (b1, k2, 0)); spec[k] += de;
      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	k = (int) sqrt (k1*k1 + k2*k2 + k3*k3);
	if (k >= K) continue;
	de = norm * (EN (k1, k2, k3) + EN (b1, k2, k3)); spec[k] += de;
	de = norm * (EN (k1, b2, k3) + EN (b1, b2, k3)); spec[k] += de;
      }
    }
  }
}
