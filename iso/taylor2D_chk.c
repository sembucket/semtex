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

int N, K, FourKon3;


int main (int    argc,
	  char** argv)
{
  CVF          U;
  CF           Work;
  Param*       Info = (Param*) calloc (1, sizeof (Param));
  complex*     Wtab;
  int          perm;
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
  printParam (stderr, Info);

  N        = Info -> ngrid;
  K        = N / 2;
  FourKon3 = (4 * K) / 3;

  /* -- Get solution from file. */

  cfield   (&U);
  readCVF  (fp, U);
  fclose   (fp);

  fprintf (stderr, "Solution energy:               %g\n", energyF (U));

  /* -- Transform to PHYSICAL space. */

  Wtab = cvector (0, K-1);
  preFFT (Wtab, K);
  for (c = 1; c <= 3; c++)
    rc3DFT (U[c], Wtab, INVERSE);

  /* -- Compute maximum velocity component. */

  fprintf (stderr, "Maximum U-velocity:            %g\n", amaxf (U[1]));
  fprintf (stderr, "Maximum V-velocity:            %g\n", amaxf (U[2]));
  fprintf (stderr, "Maximum W-velocity:            %g\n", amaxf (U[3]));

  /* -- Subtract off exact solution. */

  Taylor2D_error (U, Info, perm);

  /* -- Compute maximum error velocity component. */

  fprintf (stderr, "Maximum U-velocity error:      %g\n", amaxf (U[1]));
  fprintf (stderr, "Maximum V-velocity error:      %g\n", amaxf (U[2]));
  fprintf (stderr, "Maximum W-velocity error:      %g\n", amaxf (U[3]));

  /* -- Transform back to FOURIER space. */

  for (c = 1; c <= 3; c++) {
    rc3DFT  (U[c], Wtab, FORWARD);
    scaleFT (U[c]);
  }

  fprintf (stderr, "Error energy:                  %g\n", energyF (U));

  /* -- Output error field. */

  writeParam (stdout, Info);
  writeCVF   (stdout, U);

  return EXIT_SUCCESS;
}


