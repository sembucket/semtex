/*****************************************************************************
 * energy.c: compute solution energy.
 * 
 * $Id$
 *****************************************************************************/

#include "globals.h"

#define MAG(Z)  (Z).Re*(Z).Re + (Z).Im*(Z).Im


real  energy_Physical (const int* Dim, CVF V, const complex* Wtab)
/* ------------------------------------------------------------------------- *
 * Compute & return q^2 = <UiUi>/2: diagnostic.  Do sums in physical space.
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


real  energy (const int* Dim, const CVF V)
/* ------------------------------------------------------------------------- *
 * Compute & return q^2 = <UiUi>/2: diagnostic.  Do sums in Fourier space.
 * ------------------------------------------------------------------------- */
{
  register int       c, i;
  register real      q2;
  register complex  *u = &V[1][0][0][0],
                    *v = &V[2][0][0][0],
                    *w = &V[3][0][0][0];
  const int          Npts = Dim[1] * Dim[2] * Dim[3];

  q2 = 0.0;
  for (i = 0; i < Npts; i++) q2 += MAG (u[i]) + MAG (v[i]) + MAG (w[i]);
  q2 /= 4.0;

  return q2;
}
