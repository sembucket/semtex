/*****************************************************************************
 * nl_chk: generate velocity field for Taylor--Green vortex, operate to
 * produce nonlinear terms, subtract off analytic nonlinear terms,
 * write ISO field file to stdout.
 * Write input file details and error energy on stderr.
 *
 * usage: nl_chk  -n <cubesize> > output.fld
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

void  TaylorGreenNL_error (CVF, const int*);
void nonlinear_alt (/* input     */  CVF             U   ,
		/* output    */  CVF             G   ,
		/* workspace */  CF              F   ,
                                 CVF             WK  ,
		/* using     */  const complex*  Wtab,
                                 const complex*  Stab,
                                 const int*      Dim );

int main (int argc, char *argv[])
{
  CVF           U, G, G_old, G_temp, work;
  CF            F, F_;
  int*          Dim;
  complex*      Wtab;
  complex*      Stab;
  int           TabLen;
  Param*        U_info = (Param*) calloc (1, sizeof (Param));
  int           N, Npts, Npts_P, perm;
  FILE*         fp;
  int           cubesize;
  real          err1, err2, err3;
  register int  c, i, j, k;
  real          *u, *v, *w;

  if (argc != 3) {
    fprintf (stderr, "usage: nl_chk -n <cubesize> > output.fld\n");
    exit (EXIT_FAILURE);
  }
  
  cubesize = atoi (argv[2]);
  if (!ispow2 (cubesize)) {
    fprintf (stderr, "size must be power of 2\n");
    exit (EXIT_FAILURE);
  }

  Dim    = ivector (1, 3);
  Dim[1] = (Dim[2] = cubesize);
  Dim[3] = cubesize / 2;
  Npts   = Dim[1] * Dim[2] * Dim[3];
  Npts_P = Npts + Npts;

  allocate (&U, &G, &G_old, &work, &F, &F_, &Wtab, &Stab, Dim);
  preFFT   (Wtab, Dim[3]);
  preShift (Stab, Dim[1]);

  /* -- Compute Taylor--Green velocity field. */

  TaylorGreen (U, Dim);

  /* -- Compute maximum velocity component, then transform. */
  
  fprintf (stderr, "Maximum velocity components:        %g  %g  %g\n",
	   amaxf (U[1], Dim), amaxf (U[2], Dim), amaxf (U[3], Dim));

  for (c = 1; c <= 3; c++) {
    rc3DFT  (U[c], Dim, Wtab, FORWARD);
    scaleFT (U[c], Dim);
  }

  fprintf (stderr, "Solution energy:                    %g\n",
	   energyF (U, Dim));

  /* -- Compute nonlinear terms d(UiUj)/dxj. */

  nonlinear_alt (U, G, F, work, Wtab, Stab, Dim);

  fprintf (stderr, "Maximum nonlinear components:       %g  %g  %g\n",
	   amaxf (G[1], Dim), amaxf (G[2], Dim), amaxf (G[3], Dim));

  fprintf (stderr, "Nonlinear terms' energy:            %g\n",
	   energyF (G, Dim));

  /* -- Transform nonlinear terms to PHYSICAL space. */

  for (c = 1; c <= 3; c++) rc3DFT (G[c], Dim, Wtab, INVERSE);

  /* -- Subtract analytical solution. */

  TaylorGreenNL_error (G, Dim);

  fprintf (stderr, "Maximum nonlinear error components: %g  %g  %g\n",
	   amaxf (G[1], Dim), amaxf (G[2], Dim), amaxf (G[3], Dim));

  /* -- Transform nonlinear terms back to FOURIER space. */

  for (c = 1; c <= 3; c++) {
    rc3DFT  (G[c], Dim, Wtab, FORWARD);
    scaleFT (G[c], Dim);
  }

  fprintf (stderr, "Error energy:                       %g\n",
	   energyF (G, Dim));

  /* -- Output error field. */

  writeParam (stdout, U_info);
  writeCVF   (stdout, G, Dim);

  return EXIT_SUCCESS;
}


