/*****************************************************************************
 * div_chk.c: compute maximum divergence of an ISO file.
 *
 * Copyright (C) 1992-1999 Hugh Blackburn
 *
 * Usage: div_chk file
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

int N, K, FourKon3;


int main (int    argc,
	  char** argv)
{
  CVF      U, D;
  real**   headU;
  real**   headD;
  complex* Wtab;
  Param*   Info = (Param*) calloc (1, sizeof (Param));
  int      i, c, ntot;
  real     div, Max=0.0;
  FILE*    fp;

  if (argc != 2) {
    fprintf (stderr, "div_chk: arg count\n");
    fprintf (stderr, "Usage: div_chk file\n");
    exit    (EXIT_FAILURE);
  }

  /* -- Open file, allocate storage, read velocity Fourier components. */

  fp = efopen (argv[1], "r");

  readParam  (fp, Info);
  printParam (stdout, Info);

  N = Info -> ngrid; K = N / 2; FourKon3 = (4 * K) / 3; ntot = N * N * N;

  headU = cfield (&U);
  headD = cfield (&D);

  readCVF (fp, U);
  fclose  (fp);

  /* -- Do dUi/dXi unsummed in Fourier space. */

  deriv (U, 1, D[1], 1);
  deriv (U, 2, D[2], 2);
  deriv (U, 3, D[3], 3);

  /* -- Transform to physical space. */

  Wtab = cvector (0, K-1); 
  preFFT (Wtab, K);

  rc3DFT (D[1], Wtab, INVERSE);
  rc3DFT (D[2], Wtab, INVERSE);
  rc3DFT (D[3], Wtab, INVERSE);

  /* -- Find maximum divergence. */

  for (i = 0; i < ntot; i++) {
    div = headD[1][i] + headD[2][i] + headD[3][i];
    Max = MAX (fabs (div), Max);
  }

  printf("Maximum divergence: %g\n", Max);
  
  return (EXIT_SUCCESS);
}
