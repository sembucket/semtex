/*****************************************************************************
 * div_chk.c: compute maximum divergence of an ISO file.
 *
 * Usage: div_chk file
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"


int main (int argc, char *argv[])
{
  CVF        U, D;
  real**     headU;
  real**     headD;
  int*       Dim;
  complex*   Wtab;
  header     Run_info;
  int        i, c, Npts;
  real       div, Max=0.0;
  FILE*      fp;

  if (argc != 2) {
    fprintf(stderr, "div_chk: arg count\n");
    fprintf(stderr, "Usage: div_chk file\n");
    exit(1);
  }

  /* -- Open file, allocate storage, read velocity Fourier components. */

  fp = efopen (argv[1], "r");

  read_header (fp, &Run_info);

  Dim    = ivector (1, 3);
  Dim[1] = (Dim[2] = Run_info.N_Grid);
  Dim[3] =  Dim[1] / 2;
  Npts   = Dim[1]*Dim[2]*Dim[3];

  headU = cfield (Dim, &U);
  headD = cfield (Dim, &D);

  read_field (fp, U, Npts);
    
  fclose (fp);

  /* -- Do dUi/dXi unsummed in Fourier space. */

  deriv (U, 1, D[1], 1, Dim);
  deriv (U, 2, D[2], 2, Dim);
  deriv (U, 3, D[3], 3, Dim);

  /* -- Transform to physical space. */

  Wtab = cvector (0, Dim[3]-1); 
  preFFT (Dim[3], Wtab);

  rc3DFT (D[1], Dim, Wtab, INVERSE);
  rc3DFT (D[2], Dim, Wtab, INVERSE);
  rc3DFT (D[3], Dim, Wtab, INVERSE);

  /* -- Find maximum divergence. */

  for (i=0; i < 2*Npts; i++) {
    div = headD[1][i] + headD[2][i] + headD[3][i];
    Max = MAX (div, Max);
  }

  printf("Maximum divergence: %g\n", Max);
  
  return (EXIT_SUCCESS);
}
