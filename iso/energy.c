/*===========================================================================
 * RCS Information:
 * ----------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 *===========================================================================*/

#include "globals.h"

#define MAG(Z)  (Z).Re*(Z).Re + (Z).Im*(Z).Im


float energy_Physical(ivector Dim, complex_vector_field V, cvector Wtab)
/*===========================================================================*/
/* Compute & return q^2 = <UiUi>/2: diagnostic.  Do sums in physical space.  */
/* Also: compare with energy() below to see the discrete Parseval relation.  */
/*===========================================================================*/
{
  int    c, i, Npts;
  float  DFT_Fact, q2;
  float *u, *v, *w;
  

  u = &V[1][0][0][0].Re;    /* make fast pointers to storage */
  v = &V[2][0][0][0].Re;
  w = &V[3][0][0][0].Re;

  Npts = Dim[1]*Dim[2]*Dim[3]*2;  /* Number of points in physical space */
  
  for (c=1; c<=3; c++)            /* --> Physical space. */
    rc3DFT(V[c], Dim, Wtab, INVERSE);

  q2 = 0.0;
  for (i=0; i<Npts; i++)
       q2 += u[i]*u[i] + v[i]*v[i] + w[i]*w[i];
  q2 /= 2 * Npts;

  DFT_Fact = 2.0 / Npts;
  for (i=0; i<Npts; i++) {
    u[i] *= DFT_Fact;
    v[i] *= DFT_Fact;
    w[i] *= DFT_Fact;
  }
  for (c=1; c<=3; c++)            /* --> Fourier space. */
    rc3DFT(V[c], Dim, Wtab, FORWARD);
  
  return q2;
}



float energy(ivector Dim, complex_vector_field V)
/*===========================================================================*/
/* Compute & return q^2 = <UiUi>/2: diagnostic.  Do sums in Fourier space.   */
/*===========================================================================*/
{
  int     c, i, Npts;
  float   q2;
  complex *u, *v, *w;
  

  u = &V[1][0][0][0];             /* make fast pointers to storage */
  v = &V[2][0][0][0];
  w = &V[3][0][0][0];

  Npts = Dim[1]*Dim[2]*Dim[3];    /* Number of points in Fourier space */

  q2 = 0.0;
  for (i=0; i<Npts; i++)
       q2 += MAG(u[i]) + MAG(v[i]) + MAG(w[i]);
  q2 /= 4.0;

  return q2;
}
