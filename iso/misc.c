/*****************************************************************************
 * misc.c: miscellaneous routines.
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"


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
