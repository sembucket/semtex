/*****************************************************************************
 * energy.c: compute solution energy, enstrophy moments, other scalar
 * properties.
 * 
 * $Id$
 *****************************************************************************/

#include "iso.h"

#define MAG(Z)  (Z).Re*(Z).Re + (Z).Im*(Z).Im


real  energyP (CVF V, const complex* Wtab, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Compute & return q^2 = <UiUi>/2: diagnostic.  Do sums in PHYSICAL space.
 * Also: compare with energy() below to see the discrete Parseval relation.
 * ------------------------------------------------------------------------- */
{
  register int    i, c;
  register real   q2;
  register real  *u = &V[1][0][0][0].Re,
                 *v = &V[2][0][0][0].Re,
                 *w = &V[3][0][0][0].Re; 
  const int       Npts = 2 * Dim[1] * Dim[2] * Dim[3];
  
  for (c = 1; c <= 3; c++)            /* --> Physical space. */
    rc3DFT (V[c], Dim, Wtab, INVERSE);

  q2 = 0.0;
  for (i = 0; i < Npts; i++) q2 += u[i]*u[i] + v[i]*v[i] + w[i]*w[i];
  q2 /= 2.0 * Npts;

  for (c = 1; c <= 3; c++) {          /* --> Fourier space. */
    rc3DFT  (V[c], Dim, Wtab, FORWARD);
    scaleFT (V[c], Dim);
  }
  
  return q2;
}


real  energyF (const CVF U, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Compute & return q^2 = <UiUi>/2: diagnostic.  Do sums in FOURIER space.
 * ------------------------------------------------------------------------- */
{
  register int       i;
  register real      q2;
  register complex  *u    = &U[1][0][0][0],
                    *v    = &U[2][0][0][0],
                    *w    = &U[3][0][0][0];
  const int          Npts = Dim[1] * Dim[2] * Dim[3];

  q2 = 0.0;
  for (i = 0; i < Npts; i++) q2 += MAG (u[i]) + MAG (v[i]) + MAG (w[i]);
  q2 *= 0.25;

  return q2;
}


real  rmsEns (const CVF U, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Compute & return Omega, generalized enstrophy, of order 1.  This is just
 * the mean-squared derivative of the velocity field.  See Ref [5], eq. (3.1).
 * ------------------------------------------------------------------------- */
{
  register int       c, k1, b1, k2, b2, k3;
  register real      omega;
  register real      kSqrd;
  const    int       N        = Dim[1];
  const    int       K        = Dim[3];
  const    int       FOURKon3 = (4 * K) / 3;

  omega = 0.0;

  for (c = 1; c <= 3; c++)
    for (k1 = 1; k1 < K; k1++) {
      b1     = N - k1;
      kSqrd  = SQR (k1);
      omega += kSqrd * (MAG (U[c][k1][ 0][ 0]) +
			MAG (U[c][ 0][k1][ 0]) +
			MAG (U[c][ 0][ 0][k1]) );

      for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
	b2     = N - k2;
	kSqrd  = SQR (k1) + SQR (k2);
	omega += kSqrd * (MAG (U[c][ 0][k1][k2]) + MAG (U[c][ 0][b1][k2]) +
			  MAG (U[c][k1][ 0][k2]) + MAG (U[c][b1][ 0][k2]) +
			  MAG (U[c][k1][k2][ 0]) + MAG (U[c][b1][k2][ 0]) );

	for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	  kSqrd  = SQR (k1) + SQR (k2) + SQR (k3);
	  omega += kSqrd * (MAG (U[c][k1][k2][k3]) + MAG (U[c][b1][k2][k3]) +
			    MAG (U[c][k1][b2][k3]) + MAG (U[c][b1][b2][k3]) );

	}
      }
    }
  
  return  omega;
}


real  L2norm (const CF U, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Compute & return L2 norm = Sum U^2.  U supplied in FOURIER space.
 * ------------------------------------------------------------------------- */
{
  register int       i;
  register real      l2;
  register complex  *u    = &U[0][0][0];
  const int          Npts = Dim[1] * Dim[2] * Dim[3];

  l2 = 0.0;
  for (i = 0; i < Npts; i++) l2 += MAG (u[i]);
  l2 *= 2.0;

  return l2;
}


real  amaxf (const CF U, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Find the maximum value of scalar field U, given in PHYSICAL space.
 * ------------------------------------------------------------------------- */
{
  register int   i;
  register real  mx   = 0.0;
  register real* u    = &U[0][0][0].Re;
  const int      Npts = 2 * Dim[1] * Dim[2] * Dim[3];

  for (i = 0; i < Npts; i++) mx = MAX (fabs(u[i]), mx);

  return mx;
}
