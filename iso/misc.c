/*****************************************************************************
 * misc.c: miscellaneous routines.
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

#define PEQ(Z1,c,Z2)  (Z1)->Re += (c) * (Z2)->Re; (Z1)->Im += (c) * (Z2)->Im
#define NEQ(Z1,c,Z2)  (Z1)->Re -= (c) * (Z2)->Re; (Z1)->Im -= (c) * (Z2)->Im


void zeroVF (CVF        Z  ,
	     const int* Dim)
/* ------------------------------------------------------------------------- *
 * Zero all information in Z.
 * ------------------------------------------------------------------------- */
{
  const    int   Npts = 2 * Dim[1] * Dim[2] * Dim[3];

  memset (&Z[1][0][0][0].Re, '\0', Npts * sizeof (real));
  memset (&Z[2][0][0][0].Re, '\0', Npts * sizeof (real));
  memset (&Z[3][0][0][0].Re, '\0', Npts * sizeof (real));
}


void zeroF (CF         Z  ,
	    const int* Dim)
/* ------------------------------------------------------------------------- *
 * Zero all information in Z.
 * ------------------------------------------------------------------------- */
{
  const int Npts = 2 * Dim[1] * Dim[2] * Dim[3];

  memset (&Z[0][0][0].Re, '\0', Npts * sizeof (real));
}


void copyF (CF         tgt,
	    const CF   src,
	    const int* Dim)
/* ------------------------------------------------------------------------- *
 * Copy one scalar field to another.
 * ------------------------------------------------------------------------- */
{
  const int Npts = 2 * Dim[1] * Dim[2] * Dim[3];

  memcpy (&tgt[0][0][0].Re, &src[0][0][0].Re, Npts * sizeof (real));
}


void scaleF (CF         U    ,
	     const real alpha,
	     const int* Dim  )
/* ------------------------------------------------------------------------- *
 * Scale U by alpha.
 * ------------------------------------------------------------------------- */
{
  register int   i;
  register real* u = &U[0][0][0].Re;
  const int      Npts = 2 * Dim[1] * Dim[2] * Dim[3];

  for (i = 0; i < Npts; i++) {
    *u *= alpha;
    u++;
  }
}


void setF (CF         f1 ,
	   const CF   f2 ,
	   const int* Dim)
/* ------------------------------------------------------------------------- *
 * Set field f1 = f2, over truncated space.
 * ------------------------------------------------------------------------- */
{
  register int k1, k2, k3, b1, b2;
  const int    N = Dim[1];
  const int    K = Dim[3];
  const int    FOURKon3 = (4 * K) / 3;

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
 

void addF (CF         f1 ,
	   const CF   f2 ,
	   const int* Dim)
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
 

void subF (CF         f1 ,
	   const CF   f2 ,
	   const int* Dim)
/* ------------------------------------------------------------------------- *
 * Set field f1 -= f2, over truncated space.
 * ------------------------------------------------------------------------- */
{
  register int k1, k2, k3, b1, b2;
  const int    N = Dim[1];
  const int    K = Dim[3];
  const int    FOURKon3 = (4 * K) / 3;

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


void project (CVF        G  ,
	      CF         W  ,
	      const int* Dim)
/* ------------------------------------------------------------------------- *
 * Make divergence-free projection of G, using W as workspace.
 *
 *                       Gj = Gj - kj*(kiGi)/k^2.
 * ------------------------------------------------------------------------- */
{
  register int     k1, b1, k2, b2, k3;
  register real    c1, c2, c3, kSqrd;
  register complex *f, *g;
  const int        N        = Dim[1];
  const int        K        = Dim[3];
  const int        FOURKon3 = (4 * K) / 3;

  /* -- Use W to construct kiGi. */
  
  zeroF (W, Dim);

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;

    f = W[k1][ 0];    g = G[1][k1][ 0];    PEQ (f, k1, g);
    f = W[ 0][k1];    g = G[2][ 0][k1];    PEQ (f, k1, g);
    f = W[ 0][ 0]+k1; g = G[3][ 0][ 0]+k1; PEQ (f, k1, g);

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
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

      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {

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

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
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

      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
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
