/*****************************************************************************
 * MakeIC2: Taylor-Green vortex
 *
 * usage: makeIC2 -n <cubesize> -o <outfilename>
 *
 * Creates an initial condition file corresponding to the Taylor-Green
 * vortex.  The velocity field is given by:
 *
 *                u =  sin x cos y cos z
 *                v = -cos x sin y cos z
 *                w =  0
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


int main (int argc, char *argv[])
{
  FILE*       fp;
  char        filename[FILENAME_MAX];
  int         c, i, argnr, cubesize, Npts;
  int         paramerr = FALSE;
  CVF         IC;
  real**      head;
  int*        Dim;
  complex*    Wtab;
  Param*      Info = (Param*) calloc (1, sizeof (Param));
  real        Max_Vel;

  /* -- Process command-line arguments. */

  if (argc > 4 && argc < 8) {
    argnr = 1;
    do {
     if (argv[argnr][0] == '-') {
       if (argv[argnr][1] == 'n') {
	 argnr++;
	 cubesize = atoi(argv[argnr]);
	 if (!ispow2(cubesize)) {
	   paramerr = TRUE;
	   fprintf(stderr, "size must be power of 2\n");
	 }
       } else if (argv[argnr][1] == 'o') {
	 argnr++;
	 (void) strcpy(filename, argv[argnr]);
       } else {
	 paramerr = TRUE;
	 fprintf(stderr, "unrecognized flag: %s\n", argv[argnr]);
       }
     } else {
       paramerr = TRUE;
       fprintf(stderr, "bad flag: %s\n", argv[argnr]);
     }
     argnr++;
   } while (!(paramerr) && argnr < argc);
  } else {
    paramerr = TRUE;
    fprintf(stderr, "arg count\n");
  }
  
  if (paramerr) {
    fprintf (stderr, "Usage: makeIC -n <cubesize> -o <outfilename>\n");
    exit (EXIT_FAILURE);
  } else {
    fprintf (stderr, "Taylor-Green Initial Condition\n");
  }

  /* -- Allocate storage of IC components, zero all locations. */

  Dim    = ivector (1, 3);
  Dim[1] = (Dim[2] = cubesize);
  Dim[3] = Dim[1] / 2;
  Npts   = Dim[1] * Dim[2] * Dim[3];
  
  head = cfield  (Dim, &IC);
  Wtab = cvector (0, Dim[3]-1);

  /* -- Generate initial condition. */

  TaylorGreen (IC, Dim);

  preFFT (Wtab, Dim[3]);
  for (i = 1; i <= 3; i++) {
    rc3DFT  (IC[i], Dim, Wtab, FORWARD);
    scaleFT (IC[i], Dim);
  }

  /* -- Output to file. */

  Info -> modes   = cubesize;
  Info -> dt      = 0.01;
  Info -> step    = 0;
  Info -> Re      = 100.0;

  strcpy (Info -> name, "Taylor--Green");

  fp = efopen (filename, "w");

  writeParam (fp, Info);
  writeCVF   (fp, IC, Dim);
  fclose(fp);

  printParam  (stdout, Info, "$RCSfile$", "$Revision$");

  return EXIT_SUCCESS;
}

