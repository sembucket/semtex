/*****************************************************************************
 * misc.c: miscellaneous routines.
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

#define PEQ(Z1,c,Z2)  (Z1)->Re +=  (c) * (Z2)->Re; (Z1)->Im +=  (c) * (Z2)->Im
#define NEQ(Z1,c,Z2)  (Z1)->Re -=  (c) * (Z2)->Re; (Z1)->Im -=  (c) * (Z2)->Im


void  zeroVF (CVF Z, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Zero all information in Z.
 * ------------------------------------------------------------------------- */
{
  const    int   Npts = 2 * Dim[1] * Dim[2] * Dim[3];

  memset (&Z[1][0][0][0].Re, '\0', Npts * sizeof (real));
  memset (&Z[2][0][0][0].Re, '\0', Npts * sizeof (real));
  memset (&Z[3][0][0][0].Re, '\0', Npts * sizeof (real));
}


void  zeroF (CF Z, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Zero all information in Z.
 * ------------------------------------------------------------------------- */
{
  const    int   Npts = 2 * Dim[1] * Dim[2] * Dim[3];

  memset (&Z[0][0][0].Re, '\0', Npts * sizeof (real));
}


void  copyF (CF to, const CF from, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Copy one scalar field to another.
 * ------------------------------------------------------------------------- */
{
  const    int   Npts = 2 * Dim[1] * Dim[2] * Dim[3];

  memcpy (&to[0][0][0].Re, &from[0][0][0].Re, Npts * sizeof (real));
}


void  scaleF (CF U, const real alpha, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Scale U by alpha.
 * ------------------------------------------------------------------------- */
{
  register int    i;
  register real*  u = &U[0][0][0].Re;
  const    int    Npts = 2 * Dim[1] * Dim[2] * Dim[3];

  for (i = 0; i < Npts; i++) {
    *u *= alpha;
    u++;
  }
}


void setF (CF f1, const CF f2, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Set field f1 = f2, over truncated space.
 * ------------------------------------------------------------------------- */
{
  register int  k1, k2, k3, b1, b2;
  const    int  N = Dim[1];
  const    int  K = Dim[3];
  const    int  FOURKon3 = (4 * K) / 3;


  f1[ 0][ 0][ 0].Re = f2[ 0][ 0][ 0].Im;

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    f1[k1][ 0][ 0] = f2[k1][ 0][ 0];
    f1[ 0][k1][ 0] = f2[ 0][k1][ 0];
    f1[ 0][ 0][k1] = f2[ 0][ 0][k1];
    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2 = N - k2;
      f1[ 0][k1][k2] = f2[ 0][k1][k2];
      f1[ 0][b1][k2] = f2[ 0][b1][k2];
      f1[k1][ 0][k2] = f2[k1][ 0][k2];
      f1[b1][ 0][k2] = f2[b1][ 0][k2];
      f1[k1][k2][ 0] = f2[k1][k2][ 0];
      f1[b1][k2][ 0] = f2[b1][k2][ 0];
      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	f1[k1][k2][k3] = f2[k1][k2][k3];
	f1[b1][k2][k3] = f2[b1][k2][k3];
	f1[k1][b2][k3] = f2[k1][b2][k3];
	f1[b1][b2][k3] = f2[b1][b2][k3];
      }
    }
  }
}
 

void addF (CF f1, const CF f2, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Set field f1 += f2, over truncated space.
 * ------------------------------------------------------------------------- */
{
  register int  k1, k2, k3, b1, b2;
  const    int  N = Dim[1];
  const    int  K = Dim[3];
  const    int  FOURKon3 = (4 * K) / 3;


  f1[ 0][ 0][ 0].Re += f2[ 0][ 0][ 0].Re;

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    f1[k1][ 0][ 0].Re += f2[k1][ 0][ 0].Re;
    f1[k1][ 0][ 0].Im += f2[k1][ 0][ 0].Im;
    f1[ 0][k1][ 0].Re += f2[ 0][k1][ 0].Re;
    f1[ 0][k1][ 0].Im += f2[ 0][k1][ 0].Im;
    f1[ 0][ 0][k1].Re += f2[ 0][ 0][k1].Re;
    f1[ 0][ 0][k1].Im += f2[ 0][ 0][k1].Im;
    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
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
      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
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
 

void subF (CF f1, const CF f2, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Set field f1 -= f2, over truncated space.
 * ------------------------------------------------------------------------- */
{
  register int  k1, k2, k3, b1, b2;
  const    int  N = Dim[1];
  const    int  K = Dim[3];
  const    int  FOURKon3 = (4 * K) / 3;


  f1[ 0][ 0][ 0].Re -= f2[ 0][ 0][ 0].Re;

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    f1[k1][ 0][ 0].Re -= f2[k1][ 0][ 0].Re;
    f1[k1][ 0][ 0].Im -= f2[k1][ 0][ 0].Im;
    f1[ 0][k1][ 0].Re -= f2[ 0][k1][ 0].Re;
    f1[ 0][k1][ 0].Im -= f2[ 0][k1][ 0].Im;
    f1[ 0][ 0][k1].Re -= f2[ 0][ 0][k1].Re;
    f1[ 0][ 0][k1].Im -= f2[ 0][ 0][k1].Im;
    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
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
      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
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


void  project (CVF G,  CF W,  const int* Dim)
/* ------------------------------------------------------------------------- *
 * Make divergence-free projection of G, using W as workspace.
 *
 *                       Gj = Gj - kj*(kiGi)/k^2.
 * ------------------------------------------------------------------------- */
{
  register int  k1, b1, k2, b2, k3;
  register real kSqrd;
  const    int  N        = Dim[1];
  const    int  K        = Dim[3];
  const    int  FOURKon3 = (4 * K) / 3;

  /* -- Use W to construct kiGi. */
  
  zeroF (W, Dim);

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;

    PEQ (W[k1][ 0]+ 0, k1, G[1][k1][ 0]+ 0);
    PEQ (W[ 0][k1]+ 0, k1, G[2][ 0][k1]+ 0);
    PEQ (W[ 0][ 0]+k1, k1, G[3][ 0][ 0]+k1);

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2 = N - k2;

      PEQ (W[k1][k2]+ 0,  k1, G[1][k1][k2]+ 0); /* kk0 face. */
      PEQ (W[k1][k2]+ 0,  k2, G[2][k1][k2]+ 0);
      PEQ (W[b1][k2]+ 0, -k1, G[1][b1][k2]+ 0); /* bk0 face. */
      PEQ (W[b1][k2]+ 0,  k2, G[2][b1][k2]+ 0);
      PEQ (W[k1][ 0]+k2,  k1, G[1][k1][ 0]+k2); /* k0k face. */
      PEQ (W[k1][ 0]+k2,  k2, G[3][k1][ 0]+k2);
      PEQ (W[b1][ 0]+k2, -k1, G[1][b1][ 0]+k2); /* b0k face. */
      PEQ (W[b1][ 0]+k2,  k2, G[3][b1][ 0]+k2);
      PEQ (W[ 0][k1]+k2,  k1, G[2][ 0][k1]+k2); /* 0kk face. */
      PEQ (W[ 0][k1]+k2,  k2, G[3][ 0][k1]+k2);
      PEQ (W[ 0][b1]+k2, -k1, G[2][ 0][b1]+k2); /* 0bk face. */
      PEQ (W[ 0][b1]+k2,  k2, G[3][ 0][b1]+k2);
      
      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {

	PEQ (W[k1][k2]+k3,  k1, G[1][k1][k2]+k3);
	PEQ (W[k1][k2]+k3,  k2, G[2][k1][k2]+k3);
	PEQ (W[k1][k2]+k3,  k3, G[3][k1][k2]+k3);

	PEQ (W[b1][k2]+k3, -k1, G[1][b1][k2]+k3);
	PEQ (W[b1][k2]+k3,  k2, G[2][b1][k2]+k3);
	PEQ (W[b1][k2]+k3,  k3, G[3][b1][k2]+k3);

	PEQ (W[k1][b2]+k3,  k1, G[1][k1][b2]+k3);
	PEQ (W[k1][b2]+k3, -k2, G[2][k1][b2]+k3);
	PEQ (W[k1][b2]+k3,  k3, G[3][k1][b2]+k3);

	PEQ (W[b1][b2]+k3, -k1, G[1][b1][b2]+k3);
	PEQ (W[b1][b2]+k3, -k2, G[2][b1][b2]+k3);
	PEQ (W[b1][b2]+k3,  k3, G[3][b1][b2]+k3);
      }
    }
  }

  /* -- Make Gj = Gj - kj * W / (k*k) */

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;

    kSqrd = SQR (k1);

    NEQ (G[1][k1][ 0]+ 0, k1/kSqrd, W[k1][ 0]+ 0);
    NEQ (G[2][ 0][k1]+ 0, k1/kSqrd, W[ 0][k1]+ 0);
    NEQ (G[3][ 0][ 0]+k1, k1/kSqrd, W[ 0][ 0]+k1);

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2 = N - k2;
      kSqrd = SQR (k1) + SQR (k2);

      NEQ (G[1][k1][k2]+ 0,  k1/kSqrd, W[k1][k2]+ 0); /* kk0 face. */
      NEQ (G[2][k1][k2]+ 0,  k2/kSqrd, W[k1][k2]+ 0);
      NEQ (G[1][b1][k2]+ 0, -k1/kSqrd, W[b1][k2]+ 0); /* bk0 face. */
      NEQ (G[2][b1][k2]+ 0,  k2/kSqrd, W[b1][k2]+ 0);
      NEQ (G[1][k1][ 0]+k2,  k1/kSqrd, W[k1][ 0]+k2); /* k0k face. */
      NEQ (G[3][k1][ 0]+k2,  k2/kSqrd, W[k1][ 0]+k2);
      NEQ (G[1][b1][ 0]+k2, -k1/kSqrd, W[b1][ 0]+k2); /* b0k face. */
      NEQ (G[3][b1][ 0]+k2,  k2/kSqrd, W[b1][ 0]+k2);
      NEQ (G[2][ 0][k1]+k2,  k1/kSqrd, W[ 0][k1]+k2); /* 0kk face. */
      NEQ (G[3][ 0][k1]+k2,  k2/kSqrd, W[ 0][k1]+k2);
      NEQ (G[2][ 0][b1]+k2, -k1/kSqrd, W[ 0][b1]+k2); /* 0bk face. */
      NEQ (G[3][ 0][b1]+k2,  k2/kSqrd, W[ 0][b1]+k2);
      
      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	kSqrd = SQR (k1) + SQR (k2) + SQR (k3);

	NEQ (G[1][k1][k2]+k3,  k1/kSqrd, W[k1][k2]+k3);
	NEQ (G[2][k1][k2]+k3,  k2/kSqrd, W[k1][k2]+k3);
	NEQ (G[3][k1][k2]+k3,  k3/kSqrd, W[k1][k2]+k3);

	NEQ (G[1][b1][k2]+k3, -k1/kSqrd, W[b1][k2]+k3);
	NEQ (G[2][b1][k2]+k3,  k2/kSqrd, W[b1][k2]+k3);
	NEQ (G[3][b1][k2]+k3,  k3/kSqrd, W[b1][k2]+k3);

	NEQ (G[1][k1][b2]+k3,  k1/kSqrd, W[k1][b2]+k3);
	NEQ (G[2][k1][b2]+k3, -k2/kSqrd, W[k1][b2]+k3);
	NEQ (G[3][k1][b2]+k3,  k3/kSqrd, W[k1][b2]+k3);

	NEQ (G[1][b1][b2]+k3, -k1/kSqrd, W[b1][b2]+k3);
	NEQ (G[2][b1][b2]+k3, -k2/kSqrd, W[b1][b2]+k3);
	NEQ (G[3][b1][b2]+k3,  k3/kSqrd, W[b1][b2]+k3);
      }
    }
  }
}
