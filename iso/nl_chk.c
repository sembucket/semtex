/*****************************************************************************
 * nl_chk: generate velocity field for Taylor--Green vortex, operate to
 * produce nonlinear terms, subtract off analytic nonlinear terms,
 * write ISO field file to stdout.
 * Write input file details and error energy on stderr.
 *
 * usage: nl_chk -n <cubesize> > output.fld
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

int N, K, FourKon3;

void TaylorGreenNL_error (CVF);


int main (int    argc,
	  char** argv)
{
  CVF           U, G_temp, work;
  CVF*          G = calloc (1, sizeof (CVF*));
  CF            F, F_;
  complex*      Wtab;
  complex*      Stab;
  int           TabLen;
  Param*        Info = (Param*) calloc (1, sizeof (Param));
  int           Npts, Npts_P, perm;
  FILE*         fp;
  real          err1, err2, err3;
  register int  c, i, j, k;
  real          *u, *v, *w;

  if (argc != 3) {
    fprintf (stderr, "usage: nl_chk -n <cubesize> > output.fld\n");
    exit (EXIT_FAILURE);
  }
  
  Info -> ngrid = atoi (argv[2]);
  if (!ispow2 (Info -> ngrid)) {
    fprintf (stderr, "size must be power of 2\n");
    exit (EXIT_FAILURE);
  }

  Info -> dt     = 0.01;
  Info -> step   = 0;
  Info -> kinvis = 0.01;

  strcpy ((Info -> session =
	   malloc (sizeof (int) * strlen ("T--G") + 1)),
	  "T--G");

  /* -- Allocate storage of IC components, zero all locations. */

  N        = Info -> ngrid;
  K        = N / 2;
  FourKon3 = (4 * K) / 3;

  allocate (&U, G, 1, &work, &F, &F_, &Wtab, &Stab);
  preFFT   (Wtab, K);
  preShift (Stab, N);

  /* -- Compute Taylor--Green velocity field. */

  TaylorGreen (U);

  /* -- Compute maximum velocity component, then transform. */
  
  fprintf (stderr, "Maximum velocity components:        %g  %g  %g\n",
	   amaxF (U[1]), amaxF (U[2]), amaxF (U[3]));

  for (c = 1; c <= 3; c++) {
    rc3DFT  (U[c], Wtab, FORWARD);
    scaleFT (U[c]);
  }

  fprintf (stderr, "Solution energy:                    %g\n",
	   energyF (U));

  /* -- Compute nonlinear terms d(UiUj)/dxj. */

  nonlinear (U, G[0], F, F_, work, Wtab, Stab);


  fprintf (stderr, "Maximum nonlinear components:       %g  %g  %g\n",
	   amaxF (G[0][1]), amaxF (G[0][2]), amaxF (G[0][3]));

  fprintf (stderr, "Nonlinear terms' energy:            %g\n",
	   energyF (G[0]));

  /* -- Transform nonlinear terms to PHYSICAL space. */

  for (c = 1; c <= 3; c++) rc3DFT (G[0][c], Wtab, INVERSE);

  /* -- Subtract analytical solution. */

  TaylorGreenNL_error (G[0]);

  fprintf (stderr, "Maximum nonlinear error components: %g  %g  %g\n",
	   amaxF (G[0][1]), amaxF (G[0][2]), amaxF (G[0][3]));

  /* -- Transform nonlinear terms back to FOURIER space. */

  for (c = 1; c <= 3; c++) {
    rc3DFT  (G[0][c], Wtab, FORWARD);
    scaleFT (G[0][c]);
  }

  fprintf (stderr, "Error energy:                       %g\n",
	   energyF (G[0]));

  /* -- Output error field. */

  writeParam (stdout, Info);
  writeCVF   (stdout, G[0]);

  return EXIT_SUCCESS;
}


