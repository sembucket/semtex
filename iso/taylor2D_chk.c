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


int main (int argc, char *argv[])
{
  CVF           U;
  CF            Work;
  header        U_info;
  int*          Dim;
  complex*      Wtab;
  int           N, Npts, Npts_P, perm;
  FILE*         fp;
  real          err_max;
  register int  c, i, j, k;

  if (argc != 4) {
    fprintf (stderr,"usage: taylor2D_chk -p 0||1||2 input.fld > output.fld\n");
    exit (EXIT_FAILURE);
  }

  perm = atoi (argv[2]);
  if (perm < 0 || perm > 2) 
    message ("taylor2D_chk", "permutation code must be 0, 1 or 2", ERROR);

  fp = efopen (argv[3], "r");

  read_header  (fp,     &U_info);
  print_header (stderr, &U_info);

  /* -- Set up the problem size. */

  Dim    =  ivector (1, 3);
  Dim[1] = (N = U_info.N_Grid);
  Dim[2] =  N;
  Dim[3] =  N / 2;
  Npts   =  Dim[1] * Dim[2] * Dim[3];
  Npts_P =  Npts + Npts;

  /* -- Get solution from file. */

  cfield     (Dim, &U);
  read_field (fp, U, Npts);
  fclose     (fp);

  fprintf (stderr, "Solution energy:                  %g\n", energyF (Dim, U));

  /* -- Transform to PHYSICAL space. */

  Wtab = cvector (0, Dim[3]-1);
  preFFT (Dim[3], Wtab);
  for (c = 1; c <= 3; c++)
    rc3DFT (U[c], Dim, Wtab, INVERSE);

  /* -- Compute maximum velocity component. */
  
  err_max = 0.0;
  for (c = 1; c <= 3; c++) {
    register real *u = & U[c][0][0][0].Re;
    for (i = 0; i < Npts_P; i++)
      err_max = MAX (u[i], err_max);
  }
  fprintf (stderr, "Maximum velocity component:       %g\n", err_max);

  /* -- Subtract off exact solution. */

  Taylor2D_error (Dim, U, &U_info, perm);

  /* -- Compute maximum error velocity component. */

  err_max = 0.0;
  for (c = 1; c <= 3; c++) {
    register real *u = & U[c][0][0][0].Re;
    for (i = 0; i < Npts_P; i++)
      err_max = MAX (u[i], err_max);
  }
  fprintf (stderr, "Maximum velocity component error: %g\n", err_max);

  /* -- Transform back to FOURIER space. */

  for (c = 1; c <= 3; c++) {
    rc3DFT  (U[c], Dim, Wtab, FORWARD);
    scaleFT (U[c], Dim);
  }

  fprintf (stderr, "Error energy:                     %g\n", energyF (Dim, U));

  /* -- Output error field. */

  write_header (stdout, &U_info);
  write_field  (stdout,  U, Npts);

  return EXIT_SUCCESS;
}


