/*****************************************************************************
 * MakeIC1: generate a set of initial condition Fourier coefficients for DNS.
 *
 * Usage: makeIC1 -n <cubesize> -o <outfilename> [-s <seed>]
 *
 * Method: in this version, the Fourier coefficients are generated with a
 * "top-hat" spectrum: i.e., uniform over a range of wavenumber magnitude.
 * If K = <cubesize> / 2 (Nyquist wavenumber), then we give a uniform magni-
 * tude to E(k) over K/3 <= k <= 2K/3, where k = sqrt(k1^2 + k2^2 + k3^2).
 * The method used to generate the various components of the velocity
 * Fourier coefficients is that discussed in detail by Rogallo (1981).
 * We generate the positive k3 part of the spectrum only, since we make use
 * of the fact that the Fourier coefficients are conjugate-symmetric in the
 * simulation.
 *
 * The optional command-line argument -s <seed> can be invoked with a
 * positive integer (e.g., -s 1) to generate different random components.
 * The IC components are normalized such that q^2 = <UiUi>/2 = 1.
 * The maximum velocity component is found (in physical space), and the
 * maximum timestep allowed by the CFL condition Delta_X/Max_Vel is printed
 * on standard output.  The grid-size is calculated with the assumption that
 * the length of the sides of the computational box is 2PI.
 *
 * Files: write the data to a file (named <outfile>) which contains a header
 * structure, of which only the size of the cube is of any interest (it is
 * checked later by the simulation programme) and then the complex data for
 * the three velocity components in turn.
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include "iso.h"


int main (int argc, char *argv[])
{
  FILE*      fp;
  char       filename[STR_MAX];
  int        c, i, argnr, cubesize;
  int        paramerr = FALSE, seed = -1;
  CVF        IC;
  int*       Dim;
  complex*   Wtab;
  Param*     Info = (Param*) calloc (1, sizeof (Param));
  real       Max_Vel;


  /* -- Process command-line arguments. */

  if (argc > 4 && argc < 8) {
    argnr = 1;
    do {
      if (argv[argnr][0] == '-') {
	if (argv[argnr][1] == 'n') {
	  argnr++;
	  cubesize = atoi (argv[argnr]);
	  if (!ispow2 (cubesize)) {
	    paramerr = TRUE;
	    fprintf (stderr, "size must be power of 2\n");
	  }
	} else if (argv[argnr][1] == 'o') {
	  argnr++;
	  strcpy (filename, argv[argnr]);
	} else if (argv[argnr][1] == 's') {
	  argnr++;
	  seed = atoi (argv[argnr]);
	} else {
	  paramerr = TRUE;
	  fprintf (stderr, "unrecognized flag: %s\n", argv[argnr]);
	}
      } else {
	paramerr = TRUE;
	fprintf (stderr, "bad flag: %s\n", argv[argnr]);
      }
      argnr++;
    } while (!(paramerr) && argnr < argc);
  } else {
    paramerr = TRUE;
    fprintf (stderr, "arg count\n");
  }
  
  if (paramerr) {
    fprintf (stderr,
	     "Usage: makeIC1 -n <cubesize> -o <outfilename> [-s <seed>]\n");
    exit (EXIT_FAILURE);
  }

  /* -- Allocate storage of IC components, zero all locations. */

  Dim    = ivector (1, 3);
  Dim[1] = (Dim[2] = cubesize);
  Dim[3] = cubesize / 2;

  cfield  (Dim, &IC);
  zeroVF (IC, Dim);

  Wtab = cvector (0, Dim[3]-1);
  preFFT (Wtab, Dim[3]);  

  /* -- Generate initial condition Fourier coefficients. */

  tophat (Dim, IC, seed);

  /* -- Normalize velocities to get q^2 = 1.  */

  normalize (IC, Dim);

  /* -- Find maximum velocity component. */

  Max_Vel = 0.0;
  for (c = 1; c <= 3; c++) {
    rc3DFT (IC[c], Dim, Wtab, INVERSE);
    Max_Vel = MAX (Max_Vel, amaxf (IC[c], Dim));
    rc3DFT  (IC[c], Dim, Wtab, FORWARD);
    scaleFT (IC[c], Dim);
  }

  /* -- Output to file. */

  Info -> modes   = Dim[1];
  Info -> dt      = 0.01;
  Info -> step    = 0;
  Info -> Re      = 100.0;

  strcpy (Info -> name, "Top Hat");

  fp = efopen (filename, "w");
  writeParam  (fp, Info);
  writeCVF    (fp, IC, Dim);
  fclose      (fp);

  printParam  (stdout, Info, "$RCSfile$", "$Revision$");
  
  sprintf (filename, "        %.3e", 2.0 * M_PI / (cubesize * Max_Vel));
  message ("CFL timestep estimate", filename, REMARK);
  sprintf (filename, "            %g", energyF (IC, Dim));
  message ("Check energy: q^2",     filename, REMARK);

  return EXIT_SUCCESS;
}


  
