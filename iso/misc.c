/*****************************************************************************
 * misc.c: miscellaneous routines.
 *
 * Copyright (C) 1992, 1999 Hugh Blackburn
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

#define PEQ(Z1,c,Z2)  (Z1)->Re += (c) * (Z2)->Re; (Z1)->Im += (c) * (Z2)->Im
#define NEQ(Z1,c,Z2)  (Z1)->Re -= (c) * (Z2)->Re; (Z1)->Im -= (c) * (Z2)->Im


void zeroVF (CVF Z)
/* ------------------------------------------------------------------------- *
 * Zero all information in Z.
 * ------------------------------------------------------------------------- */
{
  const int ntot = N * N * N;

  memset (&Z[1][0][0][0].Re, '\0', ntot * sizeof (real));
  memset (&Z[2][0][0][0].Re, '\0', ntot * sizeof (real));
  memset (&Z[3][0][0][0].Re, '\0', ntot * sizeof (real));
}


void zeroF (CF Z)
/* ------------------------------------------------------------------------- *
 * Zero all information in Z.
 * ------------------------------------------------------------------------- */
{
  const int ntot = N * N * N;

  memset (&Z[0][0][0].Re, '\0', ntot * sizeof (real));
}


void copyF (CF       tgt,
	    const CF src)
/* ------------------------------------------------------------------------- *
 * Copy one scalar field to another.
 * ------------------------------------------------------------------------- */
{
  const int ntot = N * N * N;

  memcpy (&tgt[0][0][0].Re, &src[0][0][0].Re, ntot * sizeof (real));
}


void scaleF (CF         U    ,
	     const real alpha)
/* ------------------------------------------------------------------------- *
 * Scale U by alpha.
 * ------------------------------------------------------------------------- */
{
  register int   i;
  register real* u = &U[0][0][0].Re;
  const int      ntot = N * N * N;

  for (i = 0; i < ntot; i++) u[i] *= alpha;
}


void setF (CF       f1,
	   const CF f2)
/* ------------------------------------------------------------------------- *
 * Set field f1 = f2, over truncated space.
 * ------------------------------------------------------------------------- */
{
  register int k1, k2, k3, b1, b2;

  f1[ 0][ 0][ 0].Re = f2[ 0][ 0][ 0].Im;

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    f1[k1][ 0][ 0] = f2[k1][ 0][ 0];
    f1[ 0][k1][ 0] = f2[ 0][k1][ 0];
    f1[ 0][ 0][k1] = f2[ 0][ 0][k1];
    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;
      f1[ 0][k1][k2] = f2[ 0][k1][k2];
      f1[ 0][b1][k2] = f2[ 0][b1][k2];
      f1[k1][ 0][k2] = f2[k1][ 0][k2];
      f1[b1][ 0][k2] = f2[b1][ 0][k2];
      f1[k1][k2][ 0] = f2[k1][k2][ 0];
      f1[b1][k2][ 0] = f2[b1][k2][ 0];
      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	f1[k1][k2][k3] = f2[k1][k2][k3];
	f1[b1][k2][k3] = f2[b1][k2][k3];
	f1[k1][b2][k3] = f2[k1][b2][k3];
	f1[b1][b2][k3] = f2[b1][b2][k3];
      }
    }
  }
}
 

void addF (CF       f1,
	   const CF f2)
/* ------------------------------------------------------------------------- *
 * Set field f1 += f2, over truncated space.
 * ------------------------------------------------------------------------- */
{
  register int  k1, k2, k3, b1, b2;

  f1[ 0][ 0][ 0].Re += f2[ 0][ 0][ 0].Re;

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    f1[k1][ 0][ 0].Re += f2[k1][ 0][ 0].Re;
    f1[k1][ 0][ 0].Im += f2[k1][ 0][ 0].Im;
    f1[ 0][k1][ 0].Re += f2[ 0][k1][ 0].Re;
    f1[ 0][k1][ 0].Im += f2[ 0][k1][ 0].Im;
    f1[ 0][ 0][k1].Re += f2[ 0][ 0][k1].Re;
    f1[ 0][ 0][k1].Im += f2[ 0][ 0][k1].Im;
    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;
      f1[ 0][k1][k2].Re += f2[ 0][k1][k2].Re;
      f1[ 0][k1][k2].Im += f2[ 0][k1][k2].Im;
      f1[ 0][b1][k2].Re += f2[ 0][b1][k2].Re;
      f1[ 0][b1][k2].Im += f2[ 0][b1][k2].Im;
      f1[k1][ 0][k2].Re += f2[k1][ 0][k2].Re;
      f1[k1][ 0][k2].Im += f2[k1][ 0][k2].Im;
      f1[b1][ 0][k2].Re += f2[b1][ 0][k2].Re;
      f1[b1][ 0][k2].Im += f2[b1][ 0][k2].Im;
      f1[k1][k2][ 0].Re += f2[k1][k2][ 0].Re;
      f1[k1][k2][ 0].Im += f2[k1][k2][ 0].Im;
      f1[b1][k2][ 0].Re += f2[b1][k2][ 0].Re;
      f1[b1][k2][ 0].Im += f2[b1][k2][ 0].Im;
      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	f1[k1][k2][k3].Re += f2[k1][k2][k3].Re;
	f1[k1][k2][k3].Im += f2[k1][k2][k3].Im;
	f1[b1][k2][k3].Re += f2[b1][k2][k3].Re;
	f1[b1][k2][k3].Im += f2[b1][k2][k3].Im;
	f1[k1][b2][k3].Re += f2[k1][b2][k3].Re;
	f1[k1][b2][k3].Im += f2[k1][b2][k3].Im;
	f1[b1][b2][k3].Re += f2[b1][b2][k3].Re;
	f1[b1][b2][k3].Im += f2[b1][b2][k3].Im;
      }
    }
  }
}
 

void subF (CF       f1,
	   const CF f2)
/* ------------------------------------------------------------------------- *
 * Set field f1 -= f2, over truncated space.
 * ------------------------------------------------------------------------- */
{
  register int k1, k2, k3, b1, b2;

  f1[ 0][ 0][ 0].Re -= f2[ 0][ 0][ 0].Re;

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    f1[k1][ 0][ 0].Re -= f2[k1][ 0][ 0].Re;
    f1[k1][ 0][ 0].Im -= f2[k1][ 0][ 0].Im;
    f1[ 0][k1][ 0].Re -= f2[ 0][k1][ 0].Re;
    f1[ 0][k1][ 0].Im -= f2[ 0][k1][ 0].Im;
    f1[ 0][ 0][k1].Re -= f2[ 0][ 0][k1].Re;
    f1[ 0][ 0][k1].Im -= f2[ 0][ 0][k1].Im;
    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;
      f1[ 0][k1][k2].Re -= f2[ 0][k1][k2].Re;
      f1[ 0][k1][k2].Im -= f2[ 0][k1][k2].Im;
      f1[ 0][b1][k2].Re -= f2[ 0][b1][k2].Re;
      f1[ 0][b1][k2].Im -= f2[ 0][b1][k2].Im;
      f1[k1][ 0][k2].Re -= f2[k1][ 0][k2].Re;
      f1[k1][ 0][k2].Im -= f2[k1][ 0][k2].Im;
      f1[b1][ 0][k2].Re -= f2[b1][ 0][k2].Re;
      f1[b1][ 0][k2].Im -= f2[b1][ 0][k2].Im;
      f1[k1][k2][ 0].Re -= f2[k1][k2][ 0].Re;
      f1[k1][k2][ 0].Im -= f2[k1][k2][ 0].Im;
      f1[b1][k2][ 0].Re -= f2[b1][k2][ 0].Re;
      f1[b1][k2][ 0].Im -= f2[b1][k2][ 0].Im;
      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	f1[k1][k2][k3].Re -= f2[k1][k2][k3].Re;
	f1[k1][k2][k3].Im -= f2[k1][k2][k3].Im;
	f1[b1][k2][k3].Re -= f2[b1][k2][k3].Re;
	f1[b1][k2][k3].Im -= f2[b1][k2][k3].Im;
	f1[k1][b2][k3].Re -= f2[k1][b2][k3].Re;
	f1[k1][b2][k3].Im -= f2[k1][b2][k3].Im;
	f1[b1][b2][k3].Re -= f2[b1][b2][k3].Re;
	f1[b1][b2][k3].Im -= f2[b1][b2][k3].Im;
      }
    }
  }
}


void project (CVF G,
	      CF  W)
/* ------------------------------------------------------------------------- *
 * Make divergence-free projection of G, using W as workspace.
 *
 *                       Gj = Gj - kj*(kiGi)/k^2.
 * ------------------------------------------------------------------------- */
{
  register int     k1, b1, k2, b2, k3;
  register real    c1, c2, c3, kSqrd;
  register complex *f, *g;

  /* -- Use W to construct kiGi. */
  
  zeroF (W);

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;

    f = W[k1][ 0];    g = G[1][k1][ 0];    PEQ (f, k1, g);
    f = W[ 0][k1];    g = G[2][ 0][k1];    PEQ (f, k1, g);
    f = W[ 0][ 0]+k1; g = G[3][ 0][ 0]+k1; PEQ (f, k1, g);

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;

      f = W[k1][k2];    g = G[1][k1][k2];    PEQ (f,  k1, g);
      f = W[k1][k2];    g = G[2][k1][k2];    PEQ (f,  k2, g);
      f = W[b1][k2];    g = G[1][b1][k2];    PEQ (f, -k1, g);
      f = W[b1][k2];    g = G[2][b1][k2];    PEQ (f,  k2, g);
      f = W[k1][ 0]+k2; g = G[1][k1][ 0]+k2; PEQ (f,  k1, g);
      f = W[k1][ 0]+k2; g = G[3][k1][ 0]+k2; PEQ (f,  k2, g);
      f = W[b1][ 0]+k2; g = G[1][b1][ 0]+k2; PEQ (f, -k1, g);
      f = W[b1][ 0]+k2; g = G[3][b1][ 0]+k2; PEQ (f,  k2, g);
      f = W[ 0][k1]+k2; g = G[2][ 0][k1]+k2; PEQ (f,  k1, g);
      f = W[ 0][k1]+k2; g = G[3][ 0][k1]+k2; PEQ (f,  k2, g);
      f = W[ 0][b1]+k2; g = G[2][ 0][b1]+k2; PEQ (f, -k1, g);
      f = W[ 0][b1]+k2; g = G[3][ 0][b1]+k2; PEQ (f,  k2, g);

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {

	f = W[k1][k2]+k3; g = G[1][k1][k2]+k3;	PEQ (f,  k1, g);
	f = W[k1][k2]+k3; g = G[2][k1][k2]+k3;	PEQ (f,  k2, g);
	f = W[k1][k2]+k3; g = G[3][k1][k2]+k3;	PEQ (f,  k3, g);
	
	f = W[b1][k2]+k3; g = G[1][b1][k2]+k3;	PEQ (f, -k1, g);
	f = W[b1][k2]+k3; g = G[2][b1][k2]+k3;	PEQ (f,  k2, g);
	f = W[b1][k2]+k3; g = G[3][b1][k2]+k3;	PEQ (f,  k3, g);
	
	f = W[k1][b2]+k3; g = G[1][k1][b2]+k3;	PEQ (f,  k1, g);
	f = W[k1][b2]+k3; g = G[2][k1][b2]+k3;	PEQ (f, -k2, g);
	f = W[k1][b2]+k3; g = G[3][k1][b2]+k3;	PEQ (f,  k3, g);

	f = W[b1][b2]+k3; g = G[1][b1][b2]+k3;	PEQ (f, -k1, g);
	f = W[b1][b2]+k3; g = G[2][b1][b2]+k3;	PEQ (f, -k2, g);
	f = W[b1][b2]+k3; g = G[3][b1][b2]+k3;	PEQ (f,  k3, g);

      }
    }
  }

  /* -- Make Gj = Gj - kj * W / (k*k) */

  for (k1 = 1; k1 < K; k1++) {
    b1    = N - k1;
    kSqrd = SQR (k1);
    c1    = k1 / kSqrd;

    f = G[1][k1][ 0];    g = W[k1][ 0];    NEQ (f, c1, g);
    f = G[2][ 0][k1];    g = W[ 0][k1];    NEQ (f, c1, g);
    f = G[3][ 0][ 0]+k1; g = W[ 0][ 0]+k1; NEQ (f, c1, g);

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2    = N - k2;
      kSqrd = SQR (k1) + SQR (k2);
      c1    = k1/kSqrd;
      c2    = k2/kSqrd;
      c3    = -c1;

      f = G[1][k1][k2];    g = W[k1][k2];    NEQ (f, c1, g);
      f = G[2][k1][k2];    g = W[k1][k2];    NEQ (f, c2, g);
      f = G[1][b1][k2];    g = W[b1][k2];    NEQ (f, c3, g);
      f = G[2][b1][k2];    g = W[b1][k2];    NEQ (f, c2, g);
      f = G[1][k1][ 0]+k2; g = W[k1][ 0]+k2; NEQ (f, c1, g);
      f = G[3][k1][ 0]+k2; g = W[k1][ 0]+k2; NEQ (f, c2, g);
      f = G[1][b1][ 0]+k2; g = W[b1][ 0]+k2; NEQ (f, c3, g);
      f = G[3][b1][ 0]+k2; g = W[b1][ 0]+k2; NEQ (f, c2, g);
      f = G[2][ 0][k1]+k2; g = W[ 0][k1]+k2; NEQ (f, c1, g);
      f = G[3][ 0][k1]+k2; g = W[ 0][k1]+k2; NEQ (f, c2, g);
      f = G[2][ 0][b1]+k2; g = W[ 0][b1]+k2; NEQ (f, c3, g);
      f = G[3][ 0][b1]+k2; g = W[ 0][b1]+k2; NEQ (f, c2, g);

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	kSqrd = SQR (k1) + SQR (k2) + SQR (k3);
	c1    = k1/kSqrd;
	c2    = k2/kSqrd;
	c3    = k3/kSqrd;

	f = G[1][k1][k2]+k3; g = W[k1][k2]+k3; NEQ (f,  c1, g);
	f = G[2][k1][k2]+k3; g = W[k1][k2]+k3; NEQ (f,  c2, g);
	f = G[3][k1][k2]+k3; g = W[k1][k2]+k3; NEQ (f,  c3, g);

	f = G[1][b1][k2]+k3; g = W[b1][k2]+k3; NEQ (f, -c1, g);
	f = G[2][b1][k2]+k3; g = W[b1][k2]+k3; NEQ (f,  c2, g);
	f = G[3][b1][k2]+k3; g = W[b1][k2]+k3; NEQ (f,  c3, g);
	
	f = G[1][k1][b2]+k3; g = W[k1][b2]+k3; NEQ (f,  c1, g);
	f = G[2][k1][b2]+k3; g = W[k1][b2]+k3; NEQ (f, -c2, g);
	f = G[3][k1][b2]+k3; g = W[k1][b2]+k3; NEQ (f,  c3, g);

	f = G[1][b1][b2]+k3; g = W[b1][b2]+k3; NEQ (f, -c1, g);
	f = G[2][b1][b2]+k3; g = W[b1][b2]+k3; NEQ (f, -c2, g);
	f = G[3][b1][b2]+k3; g = W[b1][b2]+k3; NEQ (f,  c3, g);

      }
    }
  }
}


void filter (CF    U   ,
	     real* F_re,
	     real* F_im)
/* ------------------------------------------------------------------------- *
 * Apply a filter mask, with real & imaginary coefficients, to
 * Fourier-transformed field U.  F_re and F_im must be at least K
 * long.  The filter used in 3D space is a tensor-product of the 1D
 * complex filter F.  Routine checks to see if F is purely real (F_im
 * = 0), and carries out either real or complex multiplication as
 * appropriate.
 * ------------------------------------------------------------------------- */
{
  register int     k1, k2, k3, b1, b2;
  register real    tp, tempRe;
  register complex TP, *u = &U[0][0][0];
   
  if (F_im) {			/* -- Filter mask is complex. */
    
    TP.Re = F_re[0] * F_re[0] * F_re[0];
    TP.Im = F_im[0] * F_im[0] * F_im[0];
    SHIFT (u[rm(0,0,0)], TP);

    for (k1 = 1; k1 < K; k1++) {
      b1 = N - k1;
      TP.Re = F_re[k1] * F_re[0] * F_re[0];
      TP.Im = F_im[k1] * F_im[0] * F_im[0];
      SHIFT (u[rm(k1, 0, 0)], TP);
      SHIFT (u[rm( 0,k1, 0)], TP);
      SHIFT (u[rm( 0, 0,k1)], TP);
      for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	b2 = N - k2;
	TP.Re = F_re[k1] * F_re[k2] * F_re[0];
	TP.Im = F_im[k1] * F_im[k2] * F_im[0];
	SHIFT (u[rm( 0,k1,k2)], TP);
	SHIFT (u[rm( 0,b1,k2)], TP);
	SHIFT (u[rm(k1, 0,k2)], TP);
	SHIFT (u[rm(b1, 0,k2)], TP);
	SHIFT (u[rm(k1,k2, 0)], TP);
	SHIFT (u[rm(b1,k2, 0)], TP);
	for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	  TP.Re = F_re[k1] * F_re[k2] * F_re[k3];
	  TP.Im = F_im[k1] * F_im[k2] * F_im[k3];
	  SHIFT (u[rm(k1,k2,k3)], TP);
	  SHIFT (u[rm(b1,k2,k3)], TP);
	  SHIFT (u[rm(k1,b2,k3)], TP);
	  SHIFT (u[rm(b1,b2,k3)], TP);
	}
      }
    }

  } else {			/* -- Filter mask is real. */
    
    tp = F_re[0] * F_re[0] * F_re[0];
    U[ 0][ 0][ 0].Re *= tp;

    for (k1 = 1; k1 < K; k1++) {
      b1 = N - k1;
      tp = F_re[k1] * F_re[0] * F_re[0];
      U[k1][ 0][ 0].Re *= tp;
      U[k1][ 0][ 0].Im *= tp;
      U[ 0][k1][ 0].Re *= tp;
      U[ 0][k1][ 0].Im *= tp;
      U[ 0][ 0][k1].Re *= tp;
      U[ 0][ 0][k1].Im *= tp;
      for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	b2 = N - k2;
	tp = F_re[k1] * F_re[k2] * F_re[0];
	U[ 0][k1][k2].Re *= tp;
	U[ 0][k1][k2].Im *= tp;
	U[ 0][b1][k2].Re *= tp;
	U[ 0][b1][k2].Im *= tp;
	U[k1][ 0][k2].Re *= tp;
	U[k1][ 0][k2].Im *= tp;
	U[b1][ 0][k2].Re *= tp;
	U[b1][ 0][k2].Im *= tp;
	U[k1][k2][ 0].Re *= tp;
	U[k1][k2][ 0].Im *= tp;
	U[b1][k2][ 0].Re *= tp;
	U[b1][k2][ 0].Im *= tp;
	for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	  tp = F_re[k1] * F_re[k2] * F_re[k3];
	  U[k1][k2][k3].Re *= tp;
	  U[k1][k2][k3].Im *= tp;
	  U[b1][k2][k3].Re *= tp;
	  U[b1][k2][k3].Im *= tp;
	  U[k1][b2][k3].Re *= tp;
	  U[k1][b2][k3].Im *= tp;
	  U[b1][b2][k3].Re *= tp;
	  U[b1][b2][k3].Im *= tp;
	}
      }
    }
  }
}


void truncateF (CF Z)
/* ------------------------------------------------------------------------- *
 * Compute the isotropic truncation (filtering) of Z, according to Orszag,
 * eq. (6.7).  We do the low-frequency faces, then Nyquist planes, finally
 * all the interior data.  Method is a bit wasteful, but we only do it once.
 * Note that FOURkon3 = (4*K)/3 is always rounded down, so test >, not >=.
 * ------------------------------------------------------------------------- */
{
  register int  i, bi, j, bj, k;

  /* -- First, zero all the data in Nyquist planes and data which
   *    is from Nyquist planes & stored, packed, on k=0 face.     */

  Z[0][0][0].Im = 0.0;

  /* -- Outer parts of i, j axes. */

  for (i = N - 1; i >= K; i--) {
    Z[i][0][0].Re = (Z[i][0][0].Im = 0.0);
    Z[0][i][0].Re = (Z[0][i][0].Im = 0.0);
  }

  /* -- Upper j side of k=0 face. */

  for (i = 1; i < K; i++) {
    bi = N - i;
    for (j = 1; j < K; j++) {
      bj = N - j;
      Z[ i][bj][0].Re = (Z[ i][bj][0].Im = 0.0);
      Z[bi][bj][0].Re = (Z[bi][bj][0].Im = 0.0);
    }
  }

  /* -- Nyquist planes. */

  for (i = 0; i < N; i++)	
    for (k = 0; k < K; k++) {
      Z[i][K][k].Re = (Z[i][K][k].Im = 0.0);
      Z[K][i][k].Re = (Z[K][i][k].Im = 0.0);
    }

  /* -- Face coeffcients.  */

  for (i = 1; i < K; i++) {		/* k=0 face: account for packing. */
    bi = N - i;
    for (j = 1; j < K; j++)
      if (i+j OUTSIDE) {
	Z[ i][j][0].Re = (Z[ i][j][0].Im = 0.0);
        Z[bi][j][0].Re = (Z[bi][j][0].Im = 0.0);
      }
  }

  for (j = 1; j < K; j++) {		/* i=0 face */
    bj = N - j;
    for (k = 1; k < K; k++)
      if (j+k OUTSIDE) {
	Z[0][ j][k].Re = (Z[0][ j][k].Im = 0.0);
        Z[0][bj][k].Re = (Z[0][bj][k].Im = 0.0);
      }
  }

  for (i = 1; i < K; i++) {		/* j=0 face */
    bi = N - i;
    for (k = 1; k < K; k++)
      if (i+k OUTSIDE) {
	Z[ i][0][k].Re = (Z[ i][0][k].Im = 0.0);
        Z[bi][0][k].Re = (Z[bi][0][k].Im = 0.0);
      }
  }

  /* -- Internal. */

  for (i = 1; i < K; i++) {
    bi = N - i;
    for (j = 1; j < K; j++) {
      bj = N - j;
      for (k = 1; k < K; k++)
	if (i+j OUTSIDE || j+k OUTSIDE || i+k OUTSIDE) { 
	  Z[ i][ j][k].Re = (Z[ i][ j][k].Im = 0.0);
	  Z[ i][bj][k].Re = (Z[ i][bj][k].Im = 0.0);
	  Z[bi][ j][k].Re = (Z[bi][ j][k].Im = 0.0);
	  Z[bi][bj][k].Re = (Z[bi][bj][k].Im = 0.0);
      }
    }
  }
}


void truncateVF (CVF Z)
/* ------------------------------------------------------------------------- *
 * Truncate all components of Z
 * ------------------------------------------------------------------------- */
{
  int i;
  
  for (i = 1; i <= 3; i++) truncateF (Z[i]);
}
