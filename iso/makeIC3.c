/*****************************************************************************
 * MakeIC3: 2D Taylor flows.
 *
 * usage: makeIC3 -n <cubesize> -o <outfilename> -p <0 || 1 || 2>
 *
 * Creates an initial condition file corresponding to the Taylor flow,
 * given by:
 *   u = -cos(PI x) sin(PI y) exp(-2 PI^2 \nu t)
 *   v =  sin(PI x) cos(PI y) exp(-2 PI^2 \nu t)
 *   w =  0
 *   p = -0.25 [cos(2 PI x) + cos(2 PI y)] exp(-4 PI^2 \nu t)
 *
 * In the solution, the nonlinear and pressure terms balance out, so the
 * decay in the local components comes from the diffusion alone.  In our
 * Fourier-space method, the correct solution should be produced with or
 * without including the contribution of the incompressible projection
 * of the nonlinear terms.
 *
 * Permutation code -p selects cyclic permutation of velocity components
 * (0 <==> combination above).
 *
 * The box size is 2 pi, and the Reynolds number based on box size is
 *
 *               Re = 1 [m/s] * 1 [m] / 0.01
 *
 * Conditions are set in PHYSICAL space and transformed to FOURIER space.
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"


int main (int    argc,
	  char** argv)
{
  FILE*    fp;
  char     filename[STR_MAX];
  int      c, i, argnr, cubesize;
  int      paramerr = FALSE, seed = -1;
  int      code = 0;
  CVF      IC;
  real**   head;
  int*     Dim;
  complex* Wtab;
  Param*   Info = (Param*) calloc (1, sizeof (Param));
  real     Max_Vel;

  /* -- Process command-line arguments. */

  if (argc > 5 && argc < 9) {
    argnr = 1;
    do {
     if (argv[argnr][0] == '-') {
       if (argv[argnr][1] == 'n') {
	 argnr++;
	 cubesize = atoi (argv[argnr]);
	 if (!ispow2 (cubesize)) {
	   paramerr = TRUE;
	   fprintf(stderr, "size must be power of 2\n");
	 }
       } else if (argv[argnr][1] == 'o') {
	 argnr++;
	 (void) strcpy (filename, argv[argnr]);
       } else if (argv[argnr][1] == 'p') {
	 argnr++;
	 code = atoi (argv[argnr]);
	 if (code < 0 || code > 2) {
	   paramerr = TRUE;
	   fprintf (stderr, "permutation code must be 0, 1 or 2\n");
	 }
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
    fprintf(stderr, "arg count\n");
  }
  
  if (paramerr) {
    fprintf (stderr,
	     "Usage: makeIC3 -n <cubesize> -o <outfilename> -p <0||1||2>\n");
    exit (EXIT_FAILURE);
  } else {
    fprintf (stderr, "Taylor 2D initial condition, ");
    switch (code) {
    case 0: fprintf (stderr, "x--y, uniform in z\n"); break;
    case 1: fprintf (stderr, "y--z, uniform in x\n"); break;
    case 2: fprintf (stderr, "x--z, uniform in y\n"); break;
    default:
      break;
    }
  }

  /* -- Allocate storage of IC components, zero all locations. */

  Dim    = ivector (1, 3);
  Dim[1] = (Dim[2] = cubesize);
  Dim[3] = Dim[1] / 2;
  
  head = cfield  (Dim, &IC);
  Wtab = cvector (0, Dim[3]-1);

  /* -- Generate initial condition. */

  Taylor2D (IC, Dim, code);

  preFFT (Wtab, Dim[3]);
  for (i = 1; i <= 3; i++) {
    rc3DFT  (IC[i], Dim, Wtab, FORWARD);
    scaleFT (IC[i], Dim);
  }

  /* -- Output to file. */

  Info -> modes   = Dim[1];
  Info -> dt      = 0.01;
  Info -> step    = 0;
  Info -> Re      = 100.0;

  strcpy (Info -> name, "Taylor-2D");

  fp = efopen (filename, "w");
  writeParam  (fp, Info);
  writeCVF    (fp, IC, Dim);
  fclose      (fp);

  printParam  (stdout, Info, "$RCSfile$", "$Revision$");

  return EXIT_SUCCESS;
}


  
