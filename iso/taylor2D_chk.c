/*****************************************************************************
 * taylor2D_chk: subtract analytical solution for 2D Taylor flow from
 * ISO field file, write to stdout.
 * Write input file details and error energy on stderr.
 *
 * usage: taylor2D_chk -p 0||1||2 input.fld > output.fld
 *   (numeric value following -p selects permutation of 2D solutions). 
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"


int main (int    argc,
	  char** argv)
{
  CVF          U;
  CF           Work;
  Param*       Info = (Param*) calloc (1, sizeof (Param));
  int*         Dim;
  complex*     Wtab;
  int          N, Npts, Npts_P, perm;
  FILE*        fp;
  real         err_max;
  register int c, i, j, k;

  if (argc != 4) {
    fprintf (stderr,"usage: taylor2D_chk -p 0||1||2 input.fld > output.fld\n");
    exit (EXIT_FAILURE);
  }

  perm = atoi (argv[2]);
  if (perm < 0 || perm > 2) 
    message ("taylor2D_chk", "permutation code must be 0, 1 or 2", ERROR);

  fp = efopen (argv[3], "r");

  readParam  (fp,     Info);
  printParam (stderr, Info, "$RCSfile$","$Revision$");

  /* -- Set up the problem size. */

  Dim    =  ivector (1, 3);
  Dim[1] = (N = Info -> modes);
  Dim[2] =  N;
  Dim[3] =  N / 2;
  Npts   =  Dim[1] * Dim[2] * Dim[3];
  Npts_P =  Npts + Npts;

  /* -- Get solution from file. */

  cfield   (Dim, &U);
  readCVF  (fp, U, Dim);
  fclose   (fp);

  fprintf (stderr, "Solution energy:               %g\n", energyF (U, Dim));

  /* -- Transform to PHYSICAL space. */

  Wtab = cvector (0, Dim[3]-1);
  preFFT (Wtab, Dim[3]);
  for (c = 1; c <= 3; c++)
    rc3DFT (U[c], Dim, Wtab, INVERSE);

  /* -- Compute maximum velocity component. */

  fprintf (stderr, "Maximum U-velocity:            %g\n", amaxf (U[1], Dim));
  fprintf (stderr, "Maximum V-velocity:            %g\n", amaxf (U[2], Dim));
  fprintf (stderr, "Maximum W-velocity:            %g\n", amaxf (U[3], Dim));

  /* -- Subtract off exact solution. */

  Taylor2D_error (U, Dim, Info, perm);

  /* -- Compute maximum error velocity component. */

  fprintf (stderr, "Maximum U-velocity error:      %g\n", amaxf (U[1], Dim));
  fprintf (stderr, "Maximum V-velocity error:      %g\n", amaxf (U[2], Dim));
  fprintf (stderr, "Maximum W-velocity error:      %g\n", amaxf (U[3], Dim));

  /* -- Transform back to FOURIER space. */

  for (c = 1; c <= 3; c++) {
    rc3DFT  (U[c], Dim, Wtab, FORWARD);
    scaleFT (U[c], Dim);
  }

  fprintf (stderr, "Error energy:                  %g\n", energyF (U, Dim));

  /* -- Output error field. */

  writeParam (stdout, Info);
  writeCVF   (stdout, U, Dim);

  return EXIT_SUCCESS;
}


