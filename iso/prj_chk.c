/*****************************************************************************
 * prj_chk.c: check the projection operator.  Load a vector field with
 * random modes, truncate it, project it to a divergence-free set, and
 * then compute its divergence.
 *
 * Copyright (C) 1999 Hugh Blackburn
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

#define SIZE 64

int N, K, FourKon3;


int main (int    argc,
	  char** argv)
{
  CVF      U, D;
  real**   headU;
  real**   headD;
  complex* Wtab;
  Param*   Info = (Param*) calloc (1, sizeof (Param));
  int      i, c, ntot, seed = 1;
  real     div, Max;
  FILE*    fp;

  N = SIZE; K = N / 2; FourKon3 = (4 * K) / 3; ntot = N * N * N;

  /* -- Allocate storage. */

  headU = cfield (&U);
  headD = cfield (&D);

  printf("Field size: %1d\n", SIZE);

  /* -- Set field U to have random modes. */

  for (c = 1; c <= 3; c++) 
    for (i = 0; i < ntot; i++)
      headU[c][i] = ran2PI (&seed);

  /* -- Truncate U, compute divergence of random modes. */

  zeroVF     (D);
  truncateVF (U);

  deriv (U, 1, D[1], 1);
  deriv (U, 2, D[2], 2);
  deriv (U, 3, D[3], 3);

  Wtab = cvector (0, K-1); 
  preFFT (Wtab, K);

  rc3DFT (D[1], Wtab, INVERSE);
  rc3DFT (D[2], Wtab, INVERSE);
  rc3DFT (D[3], Wtab, INVERSE);

  for (Max = 0.0, i = 0; i < ntot; i++) {
    div = headD[1][i] + headD[2][i] + headD[3][i];
    Max = MAX (fabs (div), Max);
  }

  printf ("Peak divergence of random field : %g\n", Max);

  /* -- Carry out projection, repeat computation. */

  zeroVF    (D);		/* -- D can have rubbish in it after FFT. */
  projectVF (U, D[1]);

  deriv (U, 1, D[1], 1);
  deriv (U, 2, D[2], 2);
  deriv (U, 3, D[3], 3);

  Wtab = cvector (0, K-1); 
  preFFT (Wtab, K);

  rc3DFT (D[1], Wtab, INVERSE);
  rc3DFT (D[2], Wtab, INVERSE);
  rc3DFT (D[3], Wtab, INVERSE);

  for (Max = 0.0, i = 0; i < ntot; i++) {
    div = headD[1][i] + headD[2][i] + headD[3][i];
    Max = MAX (fabs (div), Max);
  }

  printf ("Peak divergence after projection: %g\n", Max);
  
  return (EXIT_SUCCESS);
}
