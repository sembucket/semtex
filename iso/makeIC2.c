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
  string      filename;
  int         c, i, argnr, cubesize, Npts;
  int         paramerr = FALSE, seed = -1;
  CVF         IC;
  real**      head;
  int*        Dim;
  complex*    Wtab;
  header      IC_info;
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
  Dim[3] = cubesize / 2;
  Npts   = Dim[1] * Dim[2] * Dim[3];
  
  head = cfield  (Dim, &IC);
  Wtab = cvector (0, Dim[3]-1);

  /* -- Generate initial condition. */

  TaylorGreen (Dim, IC);

  preFFT (Dim[3], Wtab);
  for (i = 1; i <= 3; i++) {
    rc3DFT  (IC[i], Dim, Wtab, FORWARD);
    scaleFT (IC[i], Dim);
  }

  /* -- Output to file. */

  IC_info.Magic    = MAGIC;
  IC_info.N_Grid   = cubesize;
  IC_info.Delta_T  = 0.01;
  IC_info.N_Save   = 0;
  IC_info.Max_Step = 0;
  IC_info.N_Step   = 0;
  IC_info.K_Visc   = 0.01;

  strcpy (IC_info.Title,"Taylor-Green");
  strcpy (IC_info.IC_File,"<nil>");

  fp = efopen(filename, "w");

  write_header (fp, &IC_info);
  write_field  (fp, IC, Npts);

  fclose(fp);

  return EXIT_SUCCESS;
}
